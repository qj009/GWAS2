KIN <- function(gencode){
  XX<-lapply(gencode,function(x)t(x))
  KK<-lapply(XX,function(x)tcrossprod(x))
  # kinship matrix
  KIN<-lapply(KK,function(x)x/mean(diag(x)))
  return(KIN)
}
