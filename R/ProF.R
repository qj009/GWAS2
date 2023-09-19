#' Allele frequency calculation function
#' @description
#' This function allows you to calculate allele frequency.
#' @details
#' Additional details...
#'
#'
#' @param GG: Genotype matrix with letter code. Rows represent SNPs, and columns represent samples. Each sample contain two columns (two alleles).

#' @returns
#' It outputs a data frame with three columns: minor allele frequency, major allele, minor allele.


#' @keywords
#' @export
#' @examples
#' ProF(GG)

# input genotype matrix: row: snp; col: sample (each sample has two columns)
ProF<-function(GG){
  nm<-dim(GG)
  # code samples
  code<-1:(0.5*nm[2])
  # set initial frequency as 0
  Prof<-matrix(0,nm[1],1)
  Profl<-matrix(NA,nm[1],2)
  # n.ad = 1 means all samples are homozyous of this snp
  n.ad<-apply(GG,1,function(x)sum(x[2*code-1]==x[2*code])/(0.5*nm[2]))
  # add two fake snps as a tester
  # inorder to
  #GG<-cbind(GG[,1],GG[,1],GG)
  #GG[,1:2]<-"A"
  #GG[2*code,2]<-"T"
  # get unique alleles for each snp
  element<-lapply(1:nm[1],function(x)unique(GG[x,]))
  # only keep tester?
  # element[-c(1:2)]<-NULL
  n.element<-lengths(element)
  # pos.0: snps with more than 2 alleles should be removed
  pos.0<-which(n.element>2 & n.ad==1)
  pos.1<-setdiff(1:nm[1],pos.0)

  if (!identical(pos.1,integer(0))){
    GG1<-as.matrix(GG[pos.1,])
    opt<-apply(GG1,1,function(x) names(which.max(table(x))))
    minor<-apply(GG1,1,function(x) names(which.min(table(x))))
    #Prof[pos.1,1]<-apply(GG1,1,function(x)sum(x==x[1])/nm[2])
    Prof[pos.1,1]<-do.call("rbind",lapply(1:length(pos.1),function(x)sum(GG1[x,]==minor[x])/nm[2]))
    Profl[pos.1,1]<- opt
    Profl[pos.1,2]<- minor
    Prof <- data.frame(Prof,Profl)
    colnames(Prof) <- c("maf", "major_allele","minor_allele")
  }
  return(Prof)
}
