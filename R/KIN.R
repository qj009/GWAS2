#' Kinship matrix calculation function
#' @description
#' This function allows you to calculate kinship matrix.
#' @details
#' Additional details...
#'
#'
#' @param gencode: Genotype matrix with numeric code. Can be calculated from function GEN.CODE(). It should be a list containing two matrix. The first matrix is additive matrix, for major allele homozygous sample is 1, 0 for heterozygous sample, -1 for minor allele homozygous sample; The second matrix is dominant matrix, 1 for heterozygous sample, 0 for homozygous sample. In each matrix, the rows represent SNPs, and the columns represent samples.

#' @returns
#' It outputs kinship matrix(s) in list format.


#' @keywords
#' @export
#' @examples
#' KIN(gencode)

KIN <- function(gencode){
  XX<-lapply(gencode,function(x)t(x))
  KK<-lapply(XX,function(x)tcrossprod(x))
  # kinship matrix
  KIN<-lapply(KK,function(x)x/mean(diag(x)))
  return(KIN)
}
