#' SNP statistical test function
#' @description
#' This function allows you to do statistical test for selected SNP with user select method. It is usually used inside the SEL.SNP() function. It depends on FIXED, MIXED and RANDOM function.
#' @details
#' Additional details...
#'

#' @param z: the input one snp genotype matrix for all samples, dim: 1*n, n is the sample counts. The rows represent samples. The columns represent SNPs.If z is NULL, then it will calculate PAR.
#' @param YFIX: Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.
#' @param KIN: Kinship matrix. It can be obtained from KIN() function.

#' @param method: Association model user want to use. It could be FIXED or RANDOM.
#' @param PAR: Initial parameters for association test.

#' @returns: the output of TEST.SCAN function depends on the user needs. If no genotype input, it will output initial parameter. If genotype input is provided, it will output statistical test results based on the method user select. The statistical test result out is a three-element list:
#' 1. wald test reuslt: wald test statistic, wald test left tail probability(log), wald test P value;
#' 2. liklihood ratio test(lrt) result: statistics, left tail probability(log), P value;
#' 3. parameters: beta, sigma2, lambda*sigma2, gamma, standand error, lrt statistics, lrt P value,wald test statistic,wald test P value.


#' @keywords cats
#' @export
#' @examples
#' cat_function()

#

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
