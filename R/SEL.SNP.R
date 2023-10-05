#' SNP association detection function
#' @description
#' This function allows you to detect target phenotype associated SNPs
#' @details
#' Additional details...
#'
#'
#' @param CP SNP map matrix containing the SNPs genomic positions . It should contain 2 columns: chromosome and base pair positions respectively.
#' @param xx genotype matrix in numeric format. Can be calculated from function GEN.CODE(). It should be a list containing two matrix (Note: the input data must be the matrix class.). The first matrix is additive matrix, for major allele homozygous sample is 1, 0 for heterozygous sample, -1 for minor allele homozygous sample; The second matrix is dominant matrix, 1 for heterozygous sample, 0 for homozygous sample. In each matrix, the rows represent SNPs, and the columns represent samples.
#' @param YFIX Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.Note: the input data must be the matrix class.
#' @param KIN Kinship matrix. It can be obtained from KIN() function.
#' @param method Association model user want to use. It could be FIXED or RANDOM.
#' @param PAR Initial parameters for association test. The default is NULL. It can be calculated through function TEST.SCAN().
#' @param Effect_SNP_number_adjust .

#' @returns
#' SEL.HAP function output a list containing five element. 1. Wald test result (chr, pos, wald test statistic, wald test left tail probability (log), wald test P value; 2.liklihood ratio test(lrt) result (chr, pos,statistics, left tail probability(log), P value)); 3. full result with 9 elements (beta, sigma2, lambda*sigma2, gamma, standand error, lrt statistics, lrt P value,wald test statistic,wald test P value); 4. independent SNPs counts; 5. adjusted P value threshold.

#' @keywords
#' @export
#' @examples
#' SEL.SNP(CP, xx, YFIX, KIN, method, PAR=NULL,Effect_SNP_number_adjust = TRUE)

#
SEL.SNP <- function(CP, xx, YFIX, KIN, method, PAR=NULL,Effect_SNP_number_adjust = TRUE){
  # WALD test and LRT test result
  r.s<-list(list(),list())
  # c.e is all results from TEST.SCAN, 9 numbers for each snps
  c.e<-list()
  XX<-lapply(xx,function(x)t(x))
  CHR<-unique(CP[,1])
  for (i in CHR){
    POS.CHR<-which(CP[,1]==i)
    for (j in POS.CHR){
      POS.S<-CP[j,2]
      XS<-lapply(t(XX),function(x)x[,j])
      XS<-do.call('cbind',XS)
      XS<-as.matrix(XS)
      # if PAR is not provided, calculate it
      if (is.null(PAR)){
        PAR<-TEST.SCAN(YFIX,NULL,KIN,method,NULL)
      }
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
  P.Value = 0.05
  if(Effect_SNP_number_adjust==TRUE){
    RT<-(C.E[,3]-C.E[,5])/C.E[,3]
    CN<-sum(RT[RT>0])
    P.Threshold<-P.Value/CN
  }else{
    CN <- nrow(xx[[1]])
    P.Threshold <- P.Value/CN
  }

  result <- list(WALD.S, LRT.S, C.E, CN, P.Threshold)
  return(result)
}

