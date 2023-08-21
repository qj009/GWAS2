# get double allele version of genotype data
# input genotype matrix: row: snps; col: sample;
# output genotype matirx: row: snp; col: sample (each sample has two columns)
GG<-function(G.P){
  #n: sample counts
  n<-ncol(G.P)
  POS<-1:n
  GG<-matrix(0,nrow(G.P),2*n)
  GG[,2*POS-1]<-G.P
  GG[,2*POS]<-G.P
  return(GG)
}
