test_that("frequency calculation works", {
   library(googledrive)
   GEN_ID <-  drive_get(as_id("1VPLHa9QWiiey4N5jaUOc516xXo22_fVi"))
   drive_download(GEN_ID, overwrite = TRUE)
   (load(file="GEN.rda"))

   GEN.GG <- GEN[,-(1:2)]
   gg<-GG(GEN.GG)
   prof <- ProF(gg)
})
