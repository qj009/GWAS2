# make fake letter code for numberical code genotype matrix
SUDO.GG <- function(gencode.a){
  gencode.a[which(gencode.a==-1)]="A"
  gencode.a[which(gencode.a==1)]="C"
  return(gencode.a)
}
