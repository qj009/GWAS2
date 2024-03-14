# Haplotype association statistical test function:internal function
# @description
# This function allows you to do statistical test for selected haplotypes.It is usually used inside the SEL.HAP() function.
# @details
# Additional details...
#

# @param HAP.X letter code geno matrix for the position defined. The rows represent samples and each sample has two rows. The columns represent SNPs.
# @param YFIX Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.Note: the input data must be the matrix class.
# @param KIN Kinship matrix. It can be obtained from KIN() function.

# @param method Association model user want to use. It could be FIXED or RANDOM.
# @param PAR Initial parameters for association test. It can be calculated through function TEST.SCAN().

# @returns SEL.HAP function output # a list with three elements:
# 1. a list with two elements:
#   1) wald test reuslt: wald test statistic, wald test left tail probability(log), wald test P value, numbers of different haplotypes detected;
#   2) liklihood ratio test(lrt) result: statistics, left tail probability(log), P value,numbers of different haplotypes detected;
# 2. z matrix: code of haplotypes;
# 3. identified haplotypes (numerical code).

#' @importFrom plyr mapvalues

test.HAP<-function(HAP.X,YFIX,KIN,method,PAR){
  n<-nrow(YFIX)
  ##if (method=="FIXED")
  ##{KK.RES<-lapply(GEN.X,function(x)tcrossprod(x))
  ##KIN<-Map("-", KK, KK.RES)
  ##}
  ##else
  ##{KIN<-KK
  ##}
  ##KIN<-list(Reduce("+",KIN))
  #KIN<-lapply(KK,function(x)x/mean(diag(x)))
  # paste selected snps together
  NAME<-apply(HAP.X,1,paste,collapse="")
  # haplotype identified
  HAP.CAT<-unique(NAME)
  # convert identified haplotype into numbers
  NAME.hap<-mapvalues(NAME, from=HAP.CAT, to=c(1:length(HAP.CAT)))

  code<-1:n
  NAME.HAP<-matrix(0,2,n)
  # haplotype of first chromatid for all sample
  NAME.HAP[1,]<-NAME.hap[2*code-1]
  # haplotype of second chromatid for all sample
  NAME.HAP[2,]<-NAME.hap[2*code]
  #prevent wrong direction concat for haplotype createion
  Min<-apply(NAME.HAP,2,min)
  Max<-apply(NAME.HAP,2,max)

  NAME.HAPLOTYPE<-paste(Min,Max,sep="")
  ELEMENT<-unique(NAME.HAPLOTYPE)
  NUM.ELEMENT<-length(ELEMENT)
  my.scan<-list(c(0,1,1),c(0,1,1))

  # test haplotype
  # z matrix: code haplotypes into 1/0 matrix, for each sample, code its haplotype as 1, other haplotype as 0
  # Fixed method: remove last haplotype column, recode the sample with last haplotype all as -1 for other haplotypes ???
    z <- NULL
  if (NUM.ELEMENT>1){
    z<-matrix(0,n,NUM.ELEMENT)
    iter<-1
    for (ii in ELEMENT){
      z[NAME.HAPLOTYPE==ii,iter]<-1
      iter<-iter+1
    }
    if (method=="FIXED"){
      z<-as.matrix(z[,-ncol(z)])
      z[NAME.HAPLOTYPE==ELEMENT[NUM.ELEMENT],]<-(-1)
    }
    my.scan<-TEST.SCAN(YFIX,z,KIN,method,PAR)
    my.scan[[3]]<-NULL
  }
  my.scan<-lapply(my.scan,function(x)c(x,NUM.ELEMENT))

  return(list(my.scan,z,ELEMENT))
}

#

