test_that("SNP-based GWAS works", {
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
   snp_scan <- SEL.SNP(gen, YFIX, KIN=kin, method, PAR=PAR)

})
