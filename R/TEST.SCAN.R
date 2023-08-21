# GWAS code, scan snps
# z is the input genotype matrix, with one snp for all samples, dim: 1*n, n is the sample counts
# if z is NULL, then it will calculate PAR
TEST.SCAN<-function(YFIX,z,KIN,method,PAR){
  if (is.null(z)){
    if (method=="FIXED"){PAR<-FIX(NULL,YFIX,KIN,NULL)}
    if (method=="RANDOM"){PAR<-MIXED(YFIX,KIN)$par}
    return(PAR)
  }else{
    if (method=="FIXED"){my.scan<-FIX(z,YFIX,KIN,PAR)}
    if (method=="RANDOM"){my.scan<-RANDOM(z,YFIX,KIN,PAR)}
    return(my.scan)
  }
}
