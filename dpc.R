############################################################################################
# Preprocessing helpers as dpclust3p cannot directly work with the supplied VCF files for now
############################################################################################
#' Create a file with the SNV loci
createLociFile = function(vcfdat, outfile, chrom_col, pos_col, ref_col, alt_col) {
    loci = vcfdat[, c(chrom_col, pos_col, ref_col, alt_col)]
    loci_file = "loci.txt"
    write.table(loci, file=outfile, sep="\t", quote=F, row.names=F, col.names=F)
}

#' Convenience function that transforms counts for alt and ref alleles into alleleCounter output
mutwt2allelecounts = function(counts.alt, counts.ref, allele.alt, allele.ref) {
output = array(0, c(length(allele.ref), 4))
  nucleotides = c("A", "C", "G", "T")
  # Propagate the alt allele counts
  nucleo.index = match(allele.alt, nucleotides)
  for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = counts.alt[i]
  }
  
  # Propagate the reference allele counts
  nucleo.index = match(allele.ref, nucleotides)
  for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = counts.ref[i]
  }
  return(output)
}

#' Dump a file with allele counts in alleleCounter output format
createAlleleCountsFile = function(vcfdat, datacol, namecol, outfile) {
    # Allelic depths for the ref and alt alleles (AD), Approximate read depth (DP)
    tumour_stat = data.frame(do.call(rbind, strsplit(as.vector(vcfdat[,datacol]), split = ":", fixed = TRUE)))
    colnames(tumour_stat) = strsplit(as.vector(unique(vcfdat[,namecol])),':')[[1]]
    # get the number of mutant reads into mutReads and the total number of reads at each mutation into totalReads, then run the next line
    # totalReads <- as.integer(as.vector(tumour_stat[,'DP']))
    mutCount =  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',2)))
    # wtCount = totalReads - mutCount
    wtCount = as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',1)))
    counts_table = mutwt2allelecounts(counts.alt=mutCount, 
                                      counts.ref=wtCount, 
                                      allele.alt=as.character(vcfdat$V5), 
                                      allele.ref=as.character(vcfdat$V4))
                                      
    output = data.frame(as.character(vcfdat[,1]), vcfdat[,2], counts_table, rowSums(counts_table))
    colnames(output) = c("#CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
                                      
    write.table(output, file=outfile, sep="\t", quote=F, row.names=F)
}


############################################################################################
# Code to wrap up a run
############################################################################################
#' Function that stores the final output in a unified format on disk
#' @param clustering A clustering result
#' @param dataset The dataset that went into clustering
#' @param most.similar.mut Vector containing for each non-sampled mutation its most similar sampled mutation. The non-sampled mutation will be assigned to the same cluster
#' @param outfiles.prefix A prefix for the filenames
#' @param assign_sampled_muts Boolean whether to assign the non-sampled mutations (Default: TRUE)
#' @param write_tree Boolean whether to write a tree to file. Not all clustering methods return a tree (Default: FALSE)
#' @author sd11
writeStandardFinalOutput = function(clustering, dataset, most.similar.mut, outfiles.prefix, subsamplenames, assign_sampled_muts=T, write_tree=F) {
  num_samples = ncol(dataset$mutCount)
  
  ########################################################################
  # Check if mutation sampling has been done, if so, unpack and assign here
  ########################################################################
  if (!is.na(most.similar.mut) && assign_sampled_muts) {
    res = unsample_mutations(dataset, clustering)
    dataset = res$dataset
    clustering = res$clustering
  }
  
  ########################################################################
  # Write out the final mutation-cluster probabilities with all mutations spiked in
  ########################################################################
  if (!is.null(clustering$all.assignment.likelihoods) & !is.na(clustering$all.assignment.likelihoods)) {
    # Fetch and drop all columns that have just zeroes
    cols_all_zero = which(sapply(1:ncol(clustering$all.assignment.likelihoods), function(i) { max(clustering$all.assignment.likelihoods[,i])==0 }))
    if (length(cols_all_zero)!=0) {
      all_assignment_likelihoods = clustering$all.assignment.likelihoods[,-cols_all_zero, drop=F]
      cluster_colnames = (1:ncol(clustering$all.assignment.likelihoods))
      print(cluster_colnames)
      cluster_colnames = (1:ncol(clustering$all.assignment.likelihoods))[-cols_all_zero]
      print(cluster_colnames)
    } else {
      all_assignment_likelihoods = clustering$all.assignment.likelihoods
      cluster_colnames = 1:ncol(clustering$all.assignment.likelihoods)
    }
    
    all_assignment_likelihoods = data.frame(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], all_assignment_likelihoods, clustering$best.node.assignments)
    colnames(all_assignment_likelihoods) = c("chr", "start", "end", paste("prob.cluster", cluster_colnames, sep="."), "most.likely.cluster")
    write.table(all_assignment_likelihoods[dataset$mutationType=="SNV",], file=paste(outfiles.prefix, "_mutationClusterLikelihoods.bed", sep=""), quote=F, row.names=F, sep="\t")
    
    if (any(dataset$mutationType=="CNA")) {
      write.table(all_assignment_likelihoods[dataset$mutationType=="CNA",], file=paste(outfiles.prefix, "_mutationClusterLikelihoodsPseudoSNV.bed", sep=""), quote=F, row.names=F, sep="\t")
    }
  }
  
  ########################################################################
  # Add the removed mutations back in
  ########################################################################
  output = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], clustering$best.node.assignments, clustering$best.assignment.likelihoods)
  save(file="temp.RData", output, dataset)
  print("Adding in removed indices")
  res = add_removed_snvs(dataset, output)
  output = res$snv_assignment_table
  mutationType = res$mutationType
  
  ########################################################################
  # Save the indices of the mutations that were not used during the analysis
  ########################################################################
  save(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""), output, clustering)
  write.table(data.frame(mut.index=dataset$removed_indices), file=paste(outfiles.prefix,"_removedMutationsIndex.txt", sep=""), row.names=F, quote=F)
  
  ########################################################################
  # Save the consensus mutation assignments
  ########################################################################
  colnames(output) = c("chr", "start", "end", "cluster", "likelihood")
  write.table(output[mutationType=="SNV",], file=paste(outfiles.prefix, "_bestConsensusAssignments.bed", sep=""), quote=F, row.names=F, sep="\t")
  
  ########################################################################
  # Save the CNA assignments separately
  ########################################################################
  if (any(mutationType=="CNA")) {
    write.table(output[mutationType=="CNA",], file=paste(outfiles.prefix, "_bestConsensusAssignmentsPseudoSNV.bed", sep=""), quote=F, row.names=F, sep="\t")
    # Assign the CNAs to clusters using their pseudoSNV representations
    cndata = assign_cnas_to_clusters(dataset$cndata, output)
    write.table(cndata, file=paste(outfiles.prefix, "_bestCNAassignments.txt", sep=""), quote=F, row.names=F, sep="\t")
    
    if (!is.null(clustering$all.assignment.likelihoods) & !is.na(clustering$all.assignment.likelihoods)) {
      cna_assignment_likelihoods = get_cnas_cluster_probs(dataset$cndata, all_assignment_likelihoods[dataset$mutationType=="CNA",], colnames(all_assignment_likelihoods))
      write.table(cna_assignment_likelihoods, file=paste(outfiles.prefix, "_cnaClusterLikelihoods.bed", sep=""), quote=F, row.names=F, sep="\t")
    }
    
    # Create a new assignment table figure with the correct information
    # This removes pseudo SNVs as the assignmentTable will add an extra column for CNAs
    cluster_locations = clustering$cluster.locations
    cluster_locations[,3] = rep(0, nrow(cluster_locations))
    mut_assignments = table(output[mutationType=="SNV","cluster"])
    for (i in 1:nrow(cluster_locations)) {
      if (as.character(cluster_locations[i,1]) %in% names(mut_assignments)) {
        cluster_locations[i,3] = mut_assignments[as.character(cluster_locations[i,1])]
      }
    }
    plotAssignmentTable(cluster_locations, paste(outfiles.prefix, "_mutation_assignments.png", sep=""), cndata=cndata, num_samples=num_samples)
  } else {
    print("Saving assignment table in")
    print(paste(outfiles.prefix, "_mutation_assignments.png", sep=""))
    plotAssignmentTable(clustering$cluster.locations, paste(outfiles.prefix, "_mutation_assignments.png", sep=""), num_samples=num_samples)
    cluster_locations = clustering$cluster.locations
  }
  
  
  ########################################################################
  # Write final cluster locations
  ########################################################################
  if (ncol(clustering$cluster.locations) > 3) {
    # nD based clustering
    write.table(cluster_locations, paste(outfiles.prefix,"_bestClusterInfo.txt",sep=""), col.names=c("cluster.no", paste(samplename,subsamplenames,sep=""), "no.of.mutations"), sep="\t", quote=F, row.names=F)
  } else {
    # 1D based
    write.table(cluster_locations, paste(outfiles.prefix,"_bestClusterInfo.txt",sep=""), col.names=c("cluster.no","location","no.of.mutations"), row.names=F, sep="\t", quote=F)
  }
  
  ########################################################################
  # If tree based analysis, also save the tree
  ########################################################################
  if (write_tree) {
    write.table(clustering$best.tree, file=paste(outfiles.prefix, "_bestConsensusTree.txt", sep=""), quote=F, row.names=F, sep="\t")
  }
}

#' Add removed mutations back into the assignment table. SNVs will be assigned to the cluster of its most similar not-removed SNV
#' @param dataset A dataset object
#' @param snv_assignment_table Data frame with the mutation assignments
#' @return The snv_assignment_table with the removed mutations added into the position they were originally
#' @author sd11
add_removed_snvs = function(dataset, snv_assignment_table) {
  mutationType = as.character(dataset$mutationType)
  if (length(dataset$removed_indices) > 0) {
    for (i in dataset$removed_indices) {
      if (i==1) {
        snv_assignment_table = rbind(c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), snv_assignment_table)
        mutationType = c("SNV", mutationType)
      } else if (i >= nrow(snv_assignment_table)) {
        snv_assignment_table = rbind(snv_assignment_table, c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA))
        mutationType = c(mutationType, "SNV")
      } else {
        snv_assignment_table = rbind(snv_assignment_table[1:(i-1),], c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), snv_assignment_table[i:nrow(snv_assignment_table),])
        mutationType = c(mutationType[1:(i-1)], "SNV", mutationType[i:length(mutationType)])
      }
    }
  }
  
  # Sort the output in the same order as the dataset
  chrpos_input = paste(dataset$chromosome.not.filtered, dataset$mut.position.not.filtered, sep="_")
  if (!is.null(dataset$cndata)) {
    chrpos_input = c(chrpos_input, paste(dataset$chromosome[dataset$mutationType=="CNA"], dataset$position[dataset$mutationType=="CNA"], sep="_"))
  }
  
  chrpos_output = paste(snv_assignment_table[,1], snv_assignment_table[,3], sep="_")
  snv_assignment_table = snv_assignment_table[match(chrpos_input, chrpos_output),]
  return(list(snv_assignment_table=snv_assignment_table, mutationType=mutationType))
}


#' Assign CNA events to clusters using their pseudoSNV representation
#' @param cndata Data frame with the CNA data
#' @param snv_assignment_table Data frame with the mutation assignments, with chromosome and position-start/end expected as first three columns
#' @return The cndata object with extra column cluster_assignment
#' @author sd11
assign_cnas_to_clusters = function(cndata, snv_assignment_table) {
  cndata$cluster_assignment = "NA"
  for (i in 1:nrow(cndata)) {
    # Fetch all pseudoSNVs that represent this CNA
    selection = snv_assignment_table[,1]==as.character(cndata$chr[i]) & snv_assignment_table[,3]==as.character(cndata$startpos[i])
    if (any(selection)) {
      # Work out which cluster has the most pseudoSNVs assigned. That will be the cluster to which the CNA event is assigned
      pseudo_snvs = snv_assignment_table[selection,,drop=F]
      assign_inventory = table(pseudo_snvs[,4])
      cndata$cluster_assignment[i] = names(assign_inventory)[which.max(assign_inventory)]
    }
  }
  return(cndata)
}

#' Use the Pseudo-SNV probabilities to obtain a probability of each CNA of each cluster
#' @param cndata The copy number data data.frame
#' @param snv_assignment_likelihoods Probabilities of the pseudo-SNVs
#' @param cluster_colnames The colnames of the output bedfile to be used to select the assignment columns
#' @return A data.frame with chr,start,end,probs_per_cluster
#' @author sd11
get_cnas_cluster_probs = function(cndata, snv_assignment_likelihoods, cluster_colnames) {
  cna_cluster_probs = matrix(NA, nrow=nrow(cndata), ncol=ncol(snv_assignment_likelihoods)-4)
  for (i in 1:nrow(cndata)) {
    # Fetch all pseudoSNVs that represent this CNA
    selection = snv_assignment_likelihoods[,1]==as.character(cndata$chr[i]) & snv_assignment_likelihoods[,3]==as.character(cndata$startpos[i])
    if (any(selection)) {
      # Work out which cluster has the most pseudoSNVs assigned. That will be the cluster to which the CNA event is assigned
      pseudo_snvs = snv_assignment_likelihoods[selection,,drop=F]
      
      cna_cluster_probs[i,] = sapply(which(grepl("prob", cluster_colnames)), function(j) {
        if (any(pseudo_snvs[,j]==0)) {
          0
        } else {
          # Combine p-values using fishers' method
          pchisq(-2 * sum(log(pseudo_snvs[,j])), df=length(pseudo_snvs[,j]), lower.tail=F)
        }
      })
    }
  }
  # Obtain most likely cluster
  # First get those CNAs for which we don't have any probabilities and set them to NA
  assignments = apply(cna_cluster_probs, 1, function(x) { all(is.na(x)) })
  assignments[assignments] = NA
  # Then assign those for which we have probabilities to the cluster with the highest probability
  assignments[!is.na(assignments)] = unlist(apply(cna_cluster_probs, 1, which.max))
  
  output = data.frame(cndata[, c("chr", "startpos", "endpos", "CNA")], cna_cluster_probs, assignments)
  colnames(output) = c("chr", "startpos", "endpos", "CNA", cluster_colnames[grepl("prob", cluster_colnames)], "most.likely.cluster")
  
  return(output)
}

#' Custom DREAM challenge function that takes a Battenberg copy number profile and guesses the cellularity
#' using clonal copy number aberrations. Normally we run Battenberg, which gives the cellularity, but the
#' challenge is contains a step where we have to guess it and where Battenberg cannot be used (as its input
#' and we don't have access to the BAMs). This function calculates a cellularity guess from each clonal
#' CNA and then takes the median.
#' @param battenberg_subclones_file String that points to a Battenberg subclones output file
#' @return A cellularity guess based on clonal CNAs
#' @author sd11
guessCellularityFromClonalCopynumber = function(battenberg_subclones_file) {
  cnprofile = read.table(battenberg_subclones_file, header=T, stringsAsFactors=F)
  is_clonal_aberration = cnprofile$frac1_A==1 & cnprofile$nMaj1_A!=cnprofile$nMin1_A & cnprofile$endpos-cnprofile$startpos > 10000000
  
  if (sum(is_clonal_aberration) > 0) {
    rho_guess = sapply(which(is_clonal_aberration), function(index) {
      (2*cnprofile$BAF[index]-1) / (2*cnprofile$BAF[index]-cnprofile$BAF[index]*(cnprofile$nMaj1_A[index]+cnprofile$nMin1_A[index])-1+cnprofile$nMaj1_A[index])
    })
    return(median(rho_guess))
  } else {
    return(0)
  }
}

#' Function that builds a SNV coassignment probability matrix from the trace
#' @param mut.assignment.type Type of SNV assignment performed (i.e. what output files to expect)
#' @param dataset A dataset object
#' @param no.iters Total iterations
#' @param no.iters.burn.in Iterations to use as burnin
#' @return A nxn matrix that contains co-assignment probabilities for each pair of SNVs, including the ones not used during clustering (all 0)
#' @author sd11
get.snv.coassignment.matrix = function(mut.assignment.type, dataset, no.iters, no.iters.burn.in) {
  if (mut.assignment.type==1) {
    load("tumour_gsdata.RData")
    coassignments = mcclust::comp.psm(GS.data$S.i[no.iters.burn.in:no.iters, ])
  } else if (mut.assignment.type==4) {
    # Option 4 already has created this matrix, load it and add removed indices
    load("tumour_coassignment_matrix.RData")
  }
  snv_index = which(dataset$mutationType=="SNV")
  coassignments_snvs = coassignments[snv_index,snv_index]
  co.clustering = add.muts.back.in(coassignments_snvs, dataset$removed_indices, def.value=0)
  return(co.clustering)
}

#' Custom DREAM challenge function that takes a Battenberg copy number profile and guesses the cellularity
#' using clonal copy number aberrations. Normally we run Battenberg, which gives the cellularity, but the
#' challenge is contains a step where we have to guess it and where Battenberg cannot be used (as its input
#' and we don't have access to the BAMs). This function calculates a cellularity guess from each clonal
#' CNA and then takes the median.
#' @param battenberg_subclones_file String that points to a Battenberg subclones output file
#' @return A cellularity guess based on clonal CNAs
#' @author sd11
guessCellularityFromClonalCopynumber = function(battenberg_subclones_file) {
  cnprofile = read.table(battenberg_subclones_file, header=T, stringsAsFactors=F)
  is_clonal_aberration = cnprofile$frac1_A==1 & cnprofile$nMaj1_A!=cnprofile$nMin1_A & cnprofile$endpos-cnprofile$startpos > 10000000

  if (sum(is_clonal_aberration) > 0) {
    rho_guess = sapply(which(is_clonal_aberration), function(index) {
      (2*cnprofile$BAF[index]-1) / (2*cnprofile$BAF[index]-cnprofile$BAF[index]*(cnprofile$nMaj1_A[index]+cnprofile$nMin1_A[index])-1+cnprofile$nMaj1_A[index])
    })
    return(median(rho_guess))
  } else {
    return(0)
  }
}

writeChallengeOutput = function(battenberg_subclones_file, no.clusters, final_clusters_table, assignments, co.clustering) {
	print("Writing challenge output files")
	print("1A")
	write.table(guessCellularityFromClonalCopynumber(battenberg_subclones_file),"subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
	print("1B")
	write.table(no.clusters,"subchallenge1B.txt",row.names=F,col.names=F,quote=F,sep="\t")
	print("1C")
	write.table(final_clusters_table,"subchallenge1C.txt",row.names=F,col.names=F,quote=F,sep="\t")
	print("2A")
	write.table(assignments,"subchallenge2A.txt",row.names=F,col.names=F,quote=F,sep="\t")
	print("2B")
	write.table(co.clustering,"subchallenge2B.txt",row.names=F,col.names=F,quote=F,sep="\t")
	print("DONE")
}

############################################################################################
# Start pipeline
############################################################################################
args = commandArgs(TRUE)
vcfdat = read.table(args[1],sep='\t',comment.char='#', stringsAsFactors=F)
datacol = as.integer(args[2]) + 10
battenberg_subclones_file = toString(args[3])
battenberg_cellularity_file = toString(args[4])
coclusterCNA = F #as.logical(args[5])
mut.assignment.type = 1 #1 #as.numeric(args[6])
sex = "male"
is.male = ifelse(sex=="male", T, F)

# Set the expected chromosomes based on the sex
if (is.male) {
  supported_chroms = as.character(c(1:22, "X", "Y"))
} else {
  supported_chroms = as.character(c(1:22, "X"))
}

# Generate a representation of the rho/psi file that preprocessing needs
samplename = "tumour"
subsamples = c()
outdir = paste0(getwd(), "/")

# General parameters
iter = 1250
burn.in = 250
namecol = 9
# mut.assignment.type = 1
conc_param = 0.01
cluster_conc = 5
most.similar.mut = NA # Not allowing downsampling for now
min.mutreads = 0
min.depth = 1
min.frac.snvs = 0.01 # Minimum fraction of total SNVs required to call a subclone real

# CN co-clustering parameters
add.conflicts = F # Make the conflicts matrix in a dataset - Flag pertains to both copy number and mut2mut phasing
cna.conflicting.events.only = F # Add only those CNAs that are conflicting
num.clonal.events.to.add = 1 # Add this many clonal CNA events to the clustering
min.cna.size = 100 # Minim size in 10kb for a CNA event to be included

# Create a temp Battenberg rho_psi file for preprocessing
cellularity = read.table(battenberg_cellularity_file, header=T, stringsAsFactors=F)$cellularity
rho_psi = data.frame(rho=c(NA, cellularity, NA), psi=rep(NA, 3), distance=rep(NA, 3), is.best=rep(NA, 3))
row.names(rho_psi) = c("ASCAT", "FRAC_GENOME", "REF_SEG")
battenberg_rho_psi_file = "temp_rho_psi.txt"
write.table(rho_psi, file=battenberg_rho_psi_file, quote=F, col.names=T, row.names=T, sep="\t")
rm(rho_psi)

#########################################################################
# Perform preprocessing
#########################################################################
#packrat::on("~/R_packrat/icgc_production_v2.0/")
library(VariantAnnotation)
library(dpclust3p)
library(DPClust)

# Create loci file
loci_file = "loci.txt"
createLociFile(vcfdat, loci_file, 1,2,4,5)

# Create allelecounts file
allelecounts_file = "alleleCounts.txt"
createAlleleCountsFile(vcfdat, datacol, namecol, allelecounts_file)

dpinput_file = "allDirichletProcessInfo.txt"
# Run preprocessing
runGetDirichletProcessInfo(loci_file=loci_file, 
                           allele_frequencies_file=allelecounts_file, 
                           cellularity_file=battenberg_rho_psi_file, 
                           subclone_file=battenberg_subclones_file, 
                           gender=sex, 
                           SNP.phase.file="NA", 
                           mut.phase.file="NA", 
                           output_file=dpinput_file)

# Load the data
dataset = load.data(c(dpinput_file), 
                    cellularity=cellularity, 
                    Chromosome="chr", 
                    position="end",
                    WT.count="WT.count", 
                    mut.count="mut.count", 
                    subclonal.CN="subclonal.CN", 
                    no.chrs.bearing.mut="no.chrs.bearing.mut", 
                    mutation.copy.number="mutation.copy.number", 
                    subclonal.fraction="subclonal.fraction", 
                    phase="phase",
                    is.male=is.male,
                    is.vcf=F,
                    ref.genome.version="hg19",
                    min.depth=min.depth,
                    min.mutreads=min.mutreads,
		    supported_chroms=supported_chroms)

if (coclusterCNA) {
  # Run preprocessing on the copy number data
  cndata_file = "cnDirichletInput.txt"
  collate_bb_subclones(samplename, battenberg_subclones_file, sex, cndata_file)
  # Load the CNdata
  cndata = load.cn.data(c(cndata_file))
  # Add the cn data into the dataset
  dataset = add.in.cn.as.snv.cluster(dataset, 
                                     cndata, 
                                     add.conflicts=add.conflicts, 
                                     conflicting.events.only=cna.conflicting.events.only, 
                                     num.clonal.events.to.add=num.clonal.events.to.add,
                                     min.cna.size=min.cna.size)
  dataset$cndata = cndata
}
save(file="dataset.RData", dataset)

#########################################################################
# DPClust
#########################################################################
clustering = DirichletProcessClustering(mutCount=dataset$mutCount, 
                                        WTCount=dataset$WTCount, 
                                        no.iters=iter, 
                                        no.iters.burn.in=burn.in, 
                                        cellularity=cellularity, 
                                        totalCopyNumber=dataset$totalCopyNumber, 
                                        mutation.copy.number=dataset$mutation.copy.number,
                                        copyNumberAdjustment=dataset$copyNumberAdjustment, 
                                        mutationTypes=dataset$mutationType,
                                        samplename=samplename, 
                                        subsamplesrun=subsamples,
                                        output_folder=outdir, 
                                        conc_param=conc_param, 
                                        cluster_conc=cluster_conc,
                                        mut.assignment.type=mut.assignment.type,
                                        most.similar.mut=most.similar.mut,
                                        min.frac.snvs.cluster=min.frac.snvs)

# Write out the output
outfiles.prefix = paste(samplename, "_", iter, "iters_", burn.in, "burnin", sep="")
writeStandardFinalOutput(clustering=clustering, 
                         dataset=dataset,
                         most.similar.mut=most.similar.mut,
                         outfiles.prefix=outfiles.prefix,
                         subsamplenames=subsamples,
                         assign_sampled_muts=assign_sampled_muts,
                         write_tree=F)

#########################################################################
# Challenge output
#########################################################################
# Read in the final output to produce the required data
# outfiles.prefix = "tumour_1250iters_250burnin"
final_clusters_table = read.table(paste0(outfiles.prefix, "_bestClusterInfo.txt"), header=T, stringsAsFactors=F)
final_assignments = read.table(paste0(outfiles.prefix, "_bestConsensusAssignments.bed"), header=T, stringsAsFactors=F)

print(final_clusters_table)
no.muts = nrow(final_assignments)

# Screen for too small clusters and remove them
if (any(final_clusters_table$no.of.mutations < (min.frac.snvs*no.muts))) {
  # TODO: Maybe only unless supported by CNAs as well?
  rowids = which(final_clusters_table$no.of.mutations < (min.frac.snvs*no.muts))
  for (rowid in rowids) {
    clusterid = final_clusters_table$cluster.no[rowid]
    snvs_assigned = final_assignments$cluster==clusterid
    snvs_assigned[is.na(snvs_assigned)] = F
    
    # Reset the mutation assignments
    final_assignments[snvs_assigned, "cluster"] = NA
    final_assignments[snvs_assigned, "likelihood"] = NA
  }
  final_clusters_table = final_clusters_table[-rowids,]
}
assignments = final_assignments$cluster

# Convert the CCF cluster locations into CP
final_clusters_table$location = final_clusters_table$location * cellularity
no.clusters = nrow(final_clusters_table)

# Build the co-clustering matrix
# co.clustering = get.snv.coassignment.matrix(mut.assignment.type, dataset, iter, burn.in)
no.muts = length(assignments)
# cellularity = max(optima)
co.clustering = array(0,c(no.muts,no.muts))
for(c in 1:no.clusters){
  indices = which(assignments==c)
  co.clustering[indices,indices] = 1
}
diag(co.clustering) = 1

# Assign the not assigned mutations to a dummy cluster
assignments[is.na(assignments)] = 0
final_clusters_table = rbind(data.frame(cluster.no=0, no.of.mutations=sum(assignments==0), location=0),final_clusters_table)

# Renumber the clusters to satisfy the evaluator
assignments_temp = assignments
for (i in 1:nrow(final_clusters_table)) {
  clusterid = final_clusters_table$cluster.no[i]
  assignments[assignments_temp==clusterid] = i
}
final_clusters_table$cluster.no = 1:nrow(final_clusters_table)

print("Writing challenge output files")
writeChallengeOutput(battenberg_subclones_file, no.clusters, final_clusters_table, assignments, co.clustering)
q(save="no")
