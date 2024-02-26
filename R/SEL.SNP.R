#' SNP association detection function
#' @description
#' This function allows you to detect target phenotype associated SNPs
#' @details
#' Additional details...
#'
#'

#' @param GEN Letter code genotype matrix with SNP genomic position information. The rows represent SNPs. The columns represent samples and each sample has two columns. The first two columns are chromosome and base pair positions.
#' @param YFIX Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.Note: the input data must be the matrix class.
#' @param KIN Kinship matrix. It can be obtained from KIN() function.
#' @param method Association model user want to use. It could be FIXED or RANDOM.
#' @param PAR Initial parameters for association test. The default is NULL. It can be calculated through function TEST.SCAN().

#' @returns
#' SEL.SNP function output a list containing three element. 1. Wald test result (chr, pos, wald test statistic, wald test left tail probability (log), wald test P value; 2.liklihood ratio test(lrt) result (chr, pos,statistics, left tail probability(log), P value)); 3. full result with 9 elements (beta, sigma2, lambda*sigma2, gamma, standand error, wald test statistic,wald test P value,lrt statistics, lrt P value);

#' @keywords SNP-based GWAS
#' @export
#' @import plyr

#
SEL.SNP <- function(GEN, YFIX, KIN, method, PAR=NULL){
  # WALD test and LRT test result
  r.s<-list(list(),list())
  # c.e is all results from TEST.SCAN, 9 numbers for each snps
  c.e<-list()
  CP <- matrix(as.numeric(GEN[,1:2]),nrow(GEN),2)
  gen <- GEN[,-(1:2)]
  CHR<-unique(CP[,1])

  for (i in CHR){
    POS.CHR<-which(CP[,1]==i)
    for (j in POS.CHR){
      POS.S<-CP[j,2]
      print(paste0("chromosome is ", i, ", position is ",POS.S))
      gs<-as.matrix(gen[j,])
      n<-nrow(YFIX)
      code<-1:n
      NAME <- NULL
      allele <- as.character(unique(gs))
      het <- paste0(allele[1],allele[2])
      hom1 <- paste0(allele[1],allele[1])
      hom2 <- paste0(allele[2],allele[2])
      genotype <- c(hom1,het,hom2)
      for(k in code){
        NAME.k <- paste0(gs[2*k-1,],gs[2*k,])
        if(!(NAME.k %in% genotype)){
          NAME.k <- het
        }
        NAME<-c(NAME,NAME.k)
      }
      # snp type identified
      SNP.CAT<-unique(NAME)
      # convert identified snp type into numbers
      NAME.snp<-mapvalues(NAME, from=SNP.CAT, to=c(1:length(SNP.CAT)))


      # NAME.SNP<-matrix(0,2,n)
      # snp of first chromatid for all sample
      # NAME.SNP[1,]<-NAME.snp[2*code-1]
      # snp of second chromatid for all sample
      # NAME.SNP[2,]<-NAME.hap[2*code]
      # prevent wrong direction concat for haplotype createion
      # Min<-apply(NAME.SNP,2,min)
      # Max<-apply(NAME.SNP,2,max)

      #NAME.SNPFINAL<-paste(Min,Max,sep="")
      ELEMENT<-unique(NAME.snp)
      NUM.ELEMENT<-length(ELEMENT)
      my.scan<-list(c(0,1,1),c(0,1,1))

      # test snp
      # z matrix: code haplotypes into 1/0 matrix, for each sample, code its haplotype as 1, other haplotype as 0
      # Fixed method: remove last haplotype column, recode the sample with last haplotype all as -1 for other haplotypes ???

      # if PAR is not provided, calculate it
      if (is.null(PAR)){
        PAR<-TEST.SCAN(YFIX,NULL,KIN,method,NULL)
      }

      z <- NULL
      if (NUM.ELEMENT>1){
        z<-matrix(0,n,NUM.ELEMENT)
        iter<-1
        for (ii in ELEMENT){
          z[NAME.snp==ii,iter]<-1
          iter<-iter+1
        }
        if (method=="FIXED"){
          z<-as.matrix(z[,-ncol(z)])
          z[NAME.snp==ELEMENT[NUM.ELEMENT],]<-(-1)
        }
        my.scan<-TEST.SCAN(YFIX,z,KIN,method,PAR)

        c.e[[j]]<-my.scan[[3]]
        rs<-lapply(my.scan,function(x)c(i,POS.S,x))
        r.s[[1]][[j]]<-rs[[1]]
        r.s[[2]][[j]]<-rs[[2]]
      }

    }
  }

  WALD.S<-do.call('rbind',r.s[[1]])
  LRT.S<-do.call('rbind',r.s[[2]])
  C.E<-do.call('rbind',c.e)
  # P.Value = 0.05
  # if(Effect_SNP_number_adjust==TRUE){
  #   RT<-(C.E[,3]-C.E[,5])/C.E[,3]
  #   CN<-sum(RT[RT>0])
  #   P.Threshold<-P.Value/CN
  # }else{
  #   CN <- nrow(xx[[1]])
  #   P.Threshold <- P.Value/CN
  # }

  result <- list(WALD.S, LRT.S, C.E)
  return(result)
}

