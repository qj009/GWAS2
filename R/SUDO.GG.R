# make fake letter code for numberical code genotype matrix
SUDO.GG <- function(gencode){
  gencode[which(gencode==-1)]="A"
  gencode[which(gencode==1)]="C"
  return(gencode)
}
