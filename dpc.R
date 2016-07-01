
mutationBurdenToMutationCopyNumber<-function(burden,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(burden))){
    mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
    mutCopyNumber[is.nan(mutCopyNumber)]=0
    return(mutCopyNumber)
}

mutationCopyNumberToMutationBurden<-function(copyNumber,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(copyNumber))){
    burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
    burden[is.nan(burden)|(burden<0.000001)]=0.000001
    burden[burden>0.999999]=0.999999
    return(burden)
}

subclone.dirichlet.gibbs <- function(C=30, y, N, totalCopyNumber=array(1,length(y)),cellularity=1, normalCopyNumber=array(2,length(y)) , no.chrs.bearing.mut = array(1,length(y)),iter=1000) {
  # y is a vector of the number of reads reporting each variant
  # N is a vector of the number of reads in total across the base in question (in the same order as Y obviously!)
  # C is the maximum number of clusters in the Dirichlet process
  # iter is the number of iterations of the Gibbs sampler

  num.muts <- length(y)
  print(paste("num.muts=",num.muts,sep=""))

  # Hyperparameters for alpha
  A <- B <- 0.01

  # Set up data formats for recording iterations
  pi.h <- matrix(NA, nrow=iter, ncol=C)
  V.h <- matrix(1, nrow=iter, ncol=C)
  S.i <- matrix(NA, nrow=iter, ncol=num.muts)
  Pr.S <- matrix(NA, nrow=num.muts, ncol=C)
  alpha <- rep(NA, iter)
  mutBurdens <- array(NA, c(iter, C,num.muts))

  mutCopyNum = mutationBurdenToMutationCopyNumber(y/N,totalCopyNumber,cellularity,normalCopyNumber) / no.chrs.bearing.mut
  lower=min(mutCopyNum)
  upper=max(mutCopyNum)
  difference=upper-lower
  lower=lower-difference/10
  upper=upper+difference/10
  # randomise starting positions of clusters
  pi.h[1,]=runif(C,lower,upper)
  for(c in 1:C){
    mutBurdens[1,c,]=mutationCopyNumberToMutationBurden(pi.h[1,c]*no.chrs.bearing.mut,totalCopyNumber,cellularity,normalCopyNumber)
  }

  V.h[1,] <- c(rep(0.5,C-1), 1)
  S.i[1,] <- c(1, rep(0,num.muts-1))
  alpha[1] <- 1
  V.h[1:iter, C] <- rep(1, iter)

  for (m in 2:iter) {
    if (m / 100 == round(m/100)) {print(m)}

    # Update cluster allocation for each individual mutation
    for (k in 1:num.muts) {
      #use log-space to avoid problems with very high counts
      Pr.S[k,1] <- log(V.h[m-1,1]) + y[k]*log(mutBurdens[m-1,1,k]) + (N[k]-y[k])*log(1-mutBurdens[m-1,1,k])
      Pr.S[k,2:C] <- sapply(2:C, function(V, obs.y, obs.N, pi, curr.k, j) {log(V[j]) + sum(log(1-V[1:(j-1)])) + obs.y[curr.k]*log(pi[j]) + (obs.N[curr.k] - obs.y[curr.k])*log(1-pi[j])}, V=V.h[m-1,], pi=mutBurdens[m-1,,k], obs.y=y, obs.N = N, curr.k=k)

      if(sum(is.na(Pr.S[k,]))>0){
        print("err1")
        print(Pr.S[k,])
        print(V.h[m-1,])
        print(pi.h[m-1,])
      }
      Pr.S[k,] = Pr.S[k,] - max(Pr.S[k,],na.rm=T)
      Pr.S[k,]=exp(Pr.S[k,])
      if(sum(is.na(Pr.S[k,]))>0){
        print("err2")
        print(Pr.S[k,])
      }
      Pr.S[k,] <- Pr.S[k,] / sum(Pr.S[k,],na.rm=T)
      if(sum(is.na(Pr.S[k,]))>0){
        print("err3")
        print(Pr.S[k,])
      }
      Pr.S[k,is.na(Pr.S[k,])] = 0
    }

    S.i[m,] <- sapply(1:num.muts, function(Pr, k) {sum(rmultinom(1,1,Pr[k,]) * (1:length(Pr[k,])))}, Pr=Pr.S)

    # Update stick-breaking weights
    V.h[m,1:(C-1)]  <- sapply(1:(C-1), function(S, curr.m, curr.alpha, h) {rbeta(1, 1+sum(S[curr.m,] == h), curr.alpha+sum(S[curr.m,] > h))}, S=S.i, curr.m=m, curr.alpha=alpha[m-1])
    V.h[m,c(V.h[m,1:(C-1)] == 1,FALSE)] <- 0.999 # Need to prevent one stick from taking all the remaining weight

    # Update location of subclones / clusters
    countsPerCopyNum=N*mutationCopyNumberToMutationBurden(1,totalCopyNumber,cellularity,normalCopyNumber)*no.chrs.bearing.mut

    #190512 randomise unused pi.h
    pi.h[m,]=runif(C,lower,upper)

    mutBurdens[m,,]=mutBurdens[m-1,,]
    for(c in unique(S.i[m,])){
      pi.h[m,c] = rgamma(1,shape=sum(y[S.i[m,]==c]),rate=sum(countsPerCopyNum[S.i[m,]==c]))
      mutBurdens[m,c,]=mutationCopyNumberToMutationBurden(pi.h[m,c]*no.chrs.bearing.mut,totalCopyNumber,cellularity,normalCopyNumber)
    }

    # Update alpha
    alpha[m] <- rgamma(1, shape=C+A-1, rate=B-sum(log(1-V.h[m,1:(C-1)])))
  }
  return(list(S.i=S.i, V.h=V.h, pi.h=pi.h, alpha=alpha, y1=y, N1=N))
}

Gibbs.subclone.density.est <- function(GS.data, pngFile, density.file = gsub(".png","density.txt",pngFile), polygon.file = gsub(".png","polygonData.txt",pngFile), density.smooth = 0.1, post.burn.in.start = 3000, post.burn.in.stop = 10000, density.from = 0, y.max=5, mutationCopyNumber = NA,no.chrs.bearing.mut = NA) {
  print(paste("density.smooth=",density.smooth,sep=""))
  png(filename=pngFile,,width=1500,height=1000)
  # GS.data is the list output from the above function
  # density.smooth is the smoothing factor used in R's density() function
  # post.burn.in.start is the number of iterations to drop from the Gibbs sampler output to allow the estimates to equilibrate on the posterior

  xlabel = "mutation copy number"
  if(is.na(mutationCopyNumber)){
    print("No mutationCopyNumber. Using mutation burden")
    y <- GS.data$y1
    N <- GS.data$N1
    mutationCopyNumber = y/N
    xlabel = "mutation burden"
  }

  if(!is.na(no.chrs.bearing.mut)){
    mutationCopyNumber = mutationCopyNumber / no.chrs.bearing.mut
    xlabel = "fraction of tumour cells"
  }

  V.h.cols <- GS.data$V.h
  pi.h.cols <- GS.data$pi.h
  wts <- matrix(NA, nrow=dim(V.h.cols)[1], ncol=dim(V.h.cols)[2])
  wts[,1] <- V.h.cols[,1]
  wts[,2] <- V.h.cols[,2] * (1-V.h.cols[,1])
  for (i in 3:dim(wts)[2]) {wts[,i] <- apply(1-V.h.cols[,1:(i-1)], MARGIN=1, FUN=prod) * V.h.cols[,i]}

  post.ints <- matrix(NA, ncol=post.burn.in.stop - post.burn.in.start + 1, nrow=512)

  x.max = ceiling(max(mutationCopyNumber)*12)/10

  xx <- density(c(pi.h.cols[post.burn.in.start-1,]), weights=c(wts[post.burn.in.start,]) / sum(c(wts[post.burn.in.start,])), adjust=density.smooth, from=density.from, to=x.max)$x


  for (i in post.burn.in.start : post.burn.in.stop) {
    post.ints[,i - post.burn.in.start + 1] <- density(c(pi.h.cols[i-1,]), weights=c(wts[i,]) / sum(c(wts[i,])), adjust=density.smooth, from=density.from, to=x.max)$y
  }

  polygon.data = c(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.975), rev(apply(post.ints, MARGIN=1, FUN=quantile, probs=0.025)))
  if(is.na(y.max)){
    y.max=ceiling(max(polygon.data)/10)*10
  }

  par(mar = c(5,6,4,1)+0.1)
  hist(mutationCopyNumber, breaks=seq(-0.1, x.max, 0.025), col="lightgrey",freq=FALSE, xlab=xlabel,main="", ylim=c(0,y.max),cex.axis=2,cex.lab=2)
  polygon(c(xx, rev(xx)), polygon.data, border="plum4", col=cm.colors(1,alpha=0.3))

  yy = apply(post.ints, MARGIN=1, FUN=quantile, probs=0.5)

  lines(xx, yy, col="plum4", lwd=3)

  dev.off()
  print(paste("highest density is at ",xx[which.max(yy)],sep=""))
  write.table(cbind(xx,yy),density.file,sep="\t",col.names=c(gsub(" ",".",xlabel),"median.density"),row.names=F,quote=F)
  write.table(polygon.data,polygon.file,sep="\t",row.names=F,quote=F)
}

getClusterAssignments <- function(GS.data, density.file, window.size = 20, burn.in = 300) {
  mutReads = GS.data$y1
  totalReads = GS.data$N1
  no.muts = length(mutReads)

  S.i = GS.data$S.i
  V.h = GS.data$V.h
  pi.h = GS.data$pi.h

  no.iters = nrow(S.i)

  density = read.table(density.file,row.names=NULL,header=T,sep="\t")

  localOptima = NULL
  peak.indices = NULL
  for(i in (1+window.size):(nrow(density)-window.size)){
    if(density$median.density[i] == max(density$median.density[(i-window.size):(i+window.size)])){
      localOptima = c(localOptima,density[i,1])
      peak.indices = c(peak.indices,i)
    }
  }

  print("localOptima")
  print(localOptima)

  #ASSIGN mutations to clusters
  no.optima = length(localOptima)
  if(no.optima>1){
    boundary = array(NA,no.optima-1)
    mutation.preferences = array(0,c(no.muts,no.optima))
    for(i in 1:(no.optima-1)){
      min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
      min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))

      #what distance along the line between a pair of optima do we have to go to reach the minimum density
      boundary[i] = (density[max(min.indices),1] + density[min(min.indices),1])/2
    }

    sampledIters = (burn.in + 1) : no.iters
    #don't use the intitial state
    sampledIters = sampledIters[sampledIters!=1]
    if(length(sampledIters) > 1000){
      sampledIters=floor(burn.in + (1:1000) * (no.iters - burn.in)/1000)
    }

    S.i = data.matrix(S.i)
    for(s in sampledIters){
      temp.preferences = array(0,c(no.muts,no.optima))
      for(c in unique(S.i[s,])){
        bestOptimum = sum(pi.h[s,c]>boundary)+1
        temp.preferences[S.i[s,]==c,bestOptimum] = temp.preferences[S.i[s,]==c,bestOptimum] + 1
      }
      iter.preferences = t(sapply(1:no.muts,function(p,i){as.integer(p[i,]==max(p[i,])) / sum(p[i,]==max(p[i,]))},p=temp.preferences))
      mutation.preferences = mutation.preferences + iter.preferences
    }
    mutation.preferences = mutation.preferences / length(sampledIters)
    most.likely.cluster = sapply(1:no.muts,function(m,i){which.max(m[i,])},m=mutation.preferences)
  }else{
    most.likely.cluster = rep(1,no.muts)
  }

  return(list(localOptima=localOptima,most.likely.cluster=most.likely.cluster))
}

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

#' Filter a dirichlet input file to remove data that can't be used or is deemed unreliable
filterDat = function(dat, is.male, min.mutreads, min.depth) {
    not.there.wt = is.na(dat$WT.count)
    not.there.mut = is.na(dat$mut.count)
    not.there.cn = is.na(dat$subclonal.CN)
    not.there.cna = is.na(dat$no.chrs.bearing.mut)
    not.cna = dat$no.chrs.bearing.mut==0
    if (is.male) {
        not.on.supported.chrom = !(dat$chr %in% as.character(1:22))
    } else {
        not.on.supported.chrom = !(dat$chr %in% c(1:22, "X"))
    }
    
    cov.mean = mean(dat$WT.count+dat$mut.count)
    cov.std = sd(dat$WT.count+dat$mut.count)
    
    too.high.coverage = (dat$WT.count+dat$mut.count) > cov.mean+6*cov.std
    not.coverage.threshold.depth = (dat$WT.count+dat$mut.count) < min.depth
    not.coverage.threshold.mutreads = dat$mut.count < min.mutreads
    
    print(paste("Removed", sum(not.there.wt),"with missing WTCount", sep=" "))
    print(paste("Removed", sum(not.there.mut),"with missing mutCount", sep=" "))
    print(paste("Removed", sum(not.there.cn),"with missing totalCopyNumber", sep=" "))
    print(paste("Removed", sum(not.there.cna),"with missing copyNumberAdjustment", sep=" "))
    print(paste("Removed", sum(not.cna),"with zero copyNumberAdjustment", sep=" "))
    print(paste("Removed", sum(not.on.supported.chrom), "on not supported genomic regions", sep=" "))
    print(paste("Removed", sum(too.high.coverage), "mutations with coverage over",cov.mean+6*cov.std, sep=" "))
    print(paste("Removed", sum(not.coverage.threshold.depth), "mutations with less than", min.depth, "read depth", sep=" "))
    print(paste("Removed", sum(not.coverage.threshold.mutreads), "mutations with less than", min.mutreads, "mutant reads", sep=" "))
    
    select = !(not.there.wt | not.there.mut | not.there.cn | not.there.cna | not.cna | not.on.supported.chrom | too.high.coverage | not.coverage.threshold.depth | not.coverage.threshold.mutreads)
    dat = dat[select,]
    return(list(data=dat, removed_indices=which(!select)))
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


#MAIN CODE BLOCK
args = commandArgs(TRUE)
vcfdat = read.table(args[1],sep='\t',comment.char='#', stringsAsFactors=F)
datacol = as.integer(args[2]) + 10
battenberg_subclones_file = toString(args[3])
battenberg_cellularity_file = toString(args[4])
sex = "male"

# Generate a representation of the rho/psi file that preprocessing needs
cellularity = read.table(battenberg_cellularity_file, header=T, stringsAsFactors=F)$cellularity
rho_psi = data.frame(rho=c(NA, cellularity, NA), psi=rep(NA, 3), distance=rep(NA, 3), is.best=rep(NA, 3))
row.names(rho_psi) = c("ASCAT", "FRAC_GENOME", "REF_SEG")
battenberg_rho_psi_file = "temp_rho_psi.txt"
write.table(rho_psi, file=battenberg_rho_psi_file, quote=F, col.names=T, row.names=T, sep="\t")
rm(rho_psi)

iter = 1000
burn.in = 300
namecol = 9
min.frac.snvs = 0.01

#########################################################################
# Perform preprocessing
#########################################################################
library(VariantAnnotation)
library(dpclust3p)

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
dat = read.table(dpinput_file, stringsAsFactors=F, header=T)
cellularity = read.table(battenberg_cellularity_file, header=T, stringsAsFactors=F)$cellularity

# Filtering
dataset = filterDat(dat, is.male=sex=="male", min.mutreads=1, min.depth=3)
dat = dataset$dat
removed_indices = dataset$removed_indices

#########################################################################
# DPClust
#########################################################################

#Gibbs sampler
GS.data.binomial = subclone.dirichlet.gibbs(y=dat$mut.count,
                                            N=dat$WT.count+dat$mut.count,
                                            totalCopyNumber=dat$subclonal.CN,
                                            cellularity=cellularity,
                                            no.chrs.bearing.mut=dat$no.chrs.bearing.mut,
                                            iter=iter)

# Density estimator
density.file = "density.txt"
Gibbs.subclone.density.est(GS.data.binomial,
                            "DirichletProcessplotBinomial.png", 
                            density.file = density.file, 
                            post.burn.in.start = burn.in, 
                            post.burn.in.stop = iter, 
                            y.max=50,
                            mutationCopyNumber=dat$mutation.copy.number,
                            no.chrs.bearing.mut=dat$no.chrs.bearing.mut)

#cluster assignment
cluster.assignment <- getClusterAssignments(GS.data.binomial, density.file = density.file, burn.in = burn.in)
save.image(file="temp.RData")
#put outputs into required format
#remove empty clusters
occupied.clusters = sort(unique(cluster.assignment$most.likely.cluster))
no.clusters = length(occupied.clusters)
optima = cluster.assignment$localOptima[occupied.clusters]
assignments = match(cluster.assignment$most.likely.cluster,occupied.clusters)
no.muts = length(assignments)

# Spike the removed mutations back in
if (length(removed_indices) > 0) {
    complete_assignments = rep(NA, (length(assignments) + length(removed_indices)))
    # Keep a counter that marks where we are in the assignments vector
    counter = 1
    for (i in 1:length(complete_assignments)) {
        if (i %in% removed_indices) {
            # This one was removed, we don't have an assignment, so leave it as is and don't update the counter
            next
        } else {
            # We have an assignment for this mutation, so update the complete_assignments and increment the counter
            complete_assignments[i] = assignments[counter]
            counter = counter + 1
        }
    }
    assignments = complete_assignments
}
    
final_clusters_table = data.frame(cbind(1:no.clusters,table(assignments),optima))
colnames(final_clusters_table) = c("cluster.no", "no.of.mutations", "location")
# Screen for too small clusters and remove them
if (any(final_clusters_table$no.of.mutations < (min.frac.snvs*no.muts))) {
  rowids = which(final_clusters_table$no.of.mutations < (min.frac.snvs*no.muts))
  for (rowid in rowids) {
    clusterid = final_clusters_table$cluster.no[rowid]
    snvs_assigned = assignments==clusterid
    snvs_assigned[is.na(snvs_assigned)] = F
    
    # Reset the mutation assignments
    assignments[snvs_assigned] = NA
  }
  final_clusters_table = final_clusters_table[-rowids,]
}

# Convert the CCF cluster locations into CP
final_clusters_table$location = final_clusters_table$location * cellularity
no.clusters = nrow(final_clusters_table)

no.muts = length(assignments)
cellularity = max(optima)
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
