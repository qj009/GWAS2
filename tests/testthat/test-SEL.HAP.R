test_that("haplotype based GWAS works", {
  library(googledrive)
   GEN_ID <-  drive_get(as_id("1VPLHa9QWiiey4N5jaUOc516xXo22_fVi"))
   drive_download(GEN_ID, overwrite = TRUE)
   (load(file="GEN.rda"))

   Y_ID <-  drive_get(as_id("1AF3XGQr-MwsR928NRLM5X6SwYx20NcUA"))
   drive_download(Y_ID, overwrite = TRUE)
   (load(file="Y.rda"))

   GEN.GG <- GEN[,-(1:2)]
   gg<-GG(GEN.GG)
   xx<-GEN.CODE(gg)
   kin<-KIN(xx)
   CP<-matrix(as.numeric(GEN[,1:2]),nrow(GEN),2)
   gen<-cbind(CP,gg)
   method<-"RANDOM"
   YFIX <- as.matrix(Y[,2:3])
   PAR<-TEST.SCAN(YFIX,NULL,KIN=kin,method,NULL)

   P.threshold=3E-06

   RR.MULTI<-SEL.HAP(MAP.S=NULL, POS.S=NULL,
   GEN=gen, YFIX=YFIX, KIN=kin, nHap=2,
   method=method, p.threshold=P.threshold, RR0=NULL,
   TEST=c(1,2),PAR=PAR)
   # initial haplotype test result
   WALD.HapInitial<-RR.MULTI[[2]][[1]]
   LRT.HapInitial<-RR.MULTI[[2]][[2]]

   # Haplotype final result
   WALD.FINAL<-RR.MULTI[[1]][[1]]
   LRT.FINAL<-RR.MULTI[[1]][[2]]

})
