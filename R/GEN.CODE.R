#' Numeric coded genotype matrix creation function
#' @description
#' This function allows you to convert letter coded genptype matrix into a numeric coded genotype matrix.

#' @details
#' Additional details...
#'
#'
#' @param GG Genotype matrix with letter code. Rows represent SNPs, and columns represent samples. Each sample contain two columns (two alleles).Note: the input data must be the matrix class.

#' @returns
#' It outputs numeric genotype matrix(s) in list format.It should be a list containing two matrix. The first matrix is additive matrix, for major allele homozygous sample is 1, 0 for heterozygous sample, -1 for minor allele homozygous sample; The second matrix is dominant matrix, 1 for heterozygous sample, 0 for homozygous sample. In each matrix, the rows represent SNPs, and the columns represent samples.


#' @keywords
#' @export
#' @examples
#' GEN.CODE(GG)

GEN.CODE<-function(GG){
  XAD<-list(list(),list())
  #nm: dimension of GG, row: snp, col: sample *2
  nm<-dim(GG)
  #nm: sample count
  n<-0.5*nm[2]
  code<-1:n
  #initial matrix: all 1, 1 row, n columns
  x<-matrix(1,1,n)
  # find major allele for each snp
  opt<-apply(GG,1,function(x) names((table(x))))

  for (i in 1:nm[1]){
    # first allele
    AL<-GG[i,2*code-1]
    # second allele
    AR<-GG[i,2*code]
    XAD[[1]][[i]]<-x*(-1)
    XAD[[2]][[i]]<-x*0
    # heterozygous snps
    DD<-(AL!=AR)
    # major allele homozygous snps
    AA<-(AL==opt[i] & AR==opt[i])
    # the second list is 1 for heterozygous sample, 0 for homozygous sample
    XAD[[2]][[i]][1,DD]<-1
    # the first list for major allele homozygous sample is 1, 0 for heterozygous sample, -1 for minor allele homozygous sample
    XAD[[1]][[i]][1,AA]<-1
    XAD[[1]][[i]][1,DD]<-0
    # set NA as NA
    NN <- which(is.na(AL) | is.na(AR))
      if(length(NN)!=0){
    XAD[[1]][[i]][1,NN] <- NA
    XAD[[2]][[i]][1,NN] <- NA
          }
  }
  XX<-lapply(XAD,function(x)do.call('rbind',x))
  storage.mode(XX[[1]])<-"integer"
  storage.mode(XX[[2]])<-"integer"
  if (length(unique(as.numeric(XX[[2]])))==1){
    if(unique(as.numeric(XX[[2]]))==0){
      XX[[2]]<-NULL
    }
  }
  rm(XAD)
  gc()
  return(XX)
}
