SEL.SNP <- function(CP, xx, YFIX, KIN, method, PAR=NULL){
  XX<-lapply(xx,function(x)t(x))
  CHR<-unique(CP[,1])
  for (i in CHR){
    POS.CHR<-which(CP[,1]==i)
    for (j in POS.CHR){
      POS.S<-CP[j,2]
      XS<-lapply(t(XX),function(x)x[,j])
      XS<-do.call('cbind',XS)
      XS<-as.matrix(XS)
      Re.S<-TEST.SCAN(YFIX,XS,KIN,method,PAR)
      c.e[[j]]<-Re.S[[3]]
      Re.S[[3]]<-NULL
      rs<-lapply(Re.S,function(x)c(i,POS.S,x))
      r.s[[1]][[j]]<-rs[[1]]
      r.s[[2]][[j]]<-rs[[2]]
    }
  }

  WALD.S<-do.call('rbind',r.s[[1]])
  LRT.S<-do.call('rbind',r.s[[2]])
  C.E<-do.call('rbind',c.e)
  #???
  RT<-(C.E[,3]-C.E[,5])/C.E[,3]
  #???
  CN<-sum(RT[RT>0])
  P.Threshold<-P.Value/CN
  P.threshold<-1/CN
  result <- list(WALD.S, LRT.S, C.E, CN, P.Threshold)
  return(result)
}

