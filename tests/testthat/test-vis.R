test_that("visulization function work", {
  library(googledrive)
  SNP_ID <-  drive_get(as_id("1L4K1bKKDVu8Z2v74nFjRAkcKs6tlASYh"))
  drive_download(SNP_ID, overwrite = TRUE)
  (load(file="snp_file.rda"))

  hapi_ID <-  drive_get(as_id("1Pgzl_ARDzQa49mOq8TY1L0vQf6A3J8GC"))
  drive_download(hapi_ID, overwrite = TRUE)
  (load(file="hapi_file.rda"))

  hap_ID <-  drive_get(as_id("1snvV-WJ6LcCAbP47NVFNfbmTXvdOy3kA/"))
  drive_download(hap_ID, overwrite = TRUE)
  (load(file="hap_file.rda"))

  T.Name<-"LD"
  p <- vis(T.Name, snp_file, hapi_file, hap_file, sig_line=3e-06, ylim=9)

})
