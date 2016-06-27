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

args = commandArgs(T)
battenberg_subclones_file = toString(args[1])

print("Writing challenge output files")
print("1A")
write.table(guessCellularityFromClonalCopynumber(battenberg_subclones_file),"subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
q(save="no")
