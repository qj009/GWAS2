#' Biallelic genotype matrix creation function
#' @description
#' This function allows you to convert single allele letter coded matrix to biallelic genotype matrix. It can only be applied for sample which are pure homozyous.

#' @details
#' Additional details...

#' @param G.P  Genotype matrix with letter code. Rows represent SNPs, and columns represent samples. Each sample contain one columns (one allele).

#' @returns
#' It outputs biallelic letter coded genotype matrix.Rows represent SNPs, and columns represent samples. Each sample contain two columns (two alleles).


#' @keywords Biallelic genotype
#' @export


GG<-function(G.P){
  #n: sample counts
  n<-ncol(G.P)
  POS<-1:n
  G.B<-matrix(0,nrow(G.P),2*n)
  G.B[,2*POS-1]<-G.P
  G.B[,2*POS]<-G.P
  return(G.B)
}
