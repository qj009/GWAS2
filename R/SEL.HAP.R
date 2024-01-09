#' Haplotype association detection function
#' @description
#' This function allows you to detect target phenotype associated haplotypes.
#' @details
#' Additional details...
#'
#'
#' @param MAP.S SNP map matrix containing the genomic positions for SNPs the user wants to scan. It should contain 2 columns: chromosome and base pair positions respectively. It can be NULL.
#' @param POS.S SNP map index (relevant to the full SNP map file) referring the genomic position information for SNPs which the user want to scan, It should be a vector of numbers. It can be NULL.
#' @param GEN Letter code genotype matrix with SNP genomic position information. The rows represent SNPs. The columns represent samples and each sample has two columns. The first two columns are chromosome and base pair positions.
#' @param YFIX Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.Note: the input data must be the matrix class.
#' @param KIN Kinship matrix. It can be obtained from KIN() function.
#' @param nHap Initial haplotype length user want to start with. The default value is 2.
#' @param method Association model user want to use. It could be FIXED or RANDOM.
#' @param p.threshold The p value threshold for initial significant haplotype selection. It should be a loose threshold to keep more seed haplotype for extension step.
#' @param RR0 Initial significant haplotype selection result. It can be NULL.
#' @param TEST The statistical test user want to use for significance test. 1: wald test; 2: lrt test, c(1,2) for both.
#' @param PAR Initial parameters for association test. The default is NULL. It can be calculated through function TEST.SCAN().

#' @returns
#' SEL.HAP function output a list containing two element. The first element is FINAL haplotype scanning and statistical test results. The second selement is initial haplotype statistical test results.
#' Final matrix has 9 columns: chromosome, position index of haplotype start SNP, position index haplotype end SNP, position of haplotype start SNP, position haplotype end SNP, statistic,  left tail probability(log), P value, numbers of different haplotypes detected;

#' @keywords Haplotype-based GWAS
#' @export
#' @examples
#'

#
SEL.HAP<-function(MAP.S=NULL, POS.S=NULL, GEN, YFIX, KIN, nHap=2, method, p.threshold, RR0=NULL, TEST, PAR=NULL){
  # snp MAP file: chr, pos
  # matrix(as.numeric(GEN[,1:2]),nrow(GEN),2)
  MAP<-matrix(as.numeric(GEN[,1:2]),nrow(GEN),2)
  # chr: snps chr the user wants to scan
  chr<-unique(c(MAP.S[,1],MAP[POS.S,1]))
  # if MAP.S and POS.S are NULL, then scan all chromosomes provided
  if (identical(chr,numeric(0))){
    chr<-unique(MAP[,1])
  }
  # chr index list (for each chr) for full MAP
  Pos.Chr<-list()
  for (i in chr){
    Pos.Chr[[i]]<-which(MAP[,1]==i)
  }
  # position index matrix list
  # For full MAP searching: each list contain (n.i-nHap+1) index with 2 columns from Pos.Chr
  # For given MAP (chr & pos) or index, extract all given position index in full MAP file
  PoS<-list()
  if (is.null(MAP.S) & is.null(POS.S)){
    for (i in chr){
      n.i<-length(Pos.Chr[[i]])
      PoS[[i]]<-matrix(Pos.Chr[[i]][1:(n.i-nHap+1)],(n.i-nHap+1),2)
    }
  }else{
    if (!is.null(MAP.S)){
      for (i in chr){
        POS.i<-which(MAP.S[,1]==i)
        n.i<-length(POS.i)
        PoS[[i]]<-matrix(0,n.i,2)
        for (j in 1:n.i){
          PoS[[i]][j,]<-Pos.Chr[[i]][which(MAP[Pos.Chr[[i]],2] %in% range(MAP.S[POS.i[j],-1]))]
        }
      }
    }
    if (!is.null(POS.S)){
      POS.S.CHR<-MAP[POS.S[,1],1]
      for (i in chr){
        POS.i<-which(POS.S.CHR==i)
        n.i<-length(POS.i)
        PoS[[i]]<-matrix(POS.S[POS.i,],n.i,2)
      }
    }
  }

  # PoS is a list with n chromsome, each chromsome is a matirx with two columns:
  # the first column is original position index, the second column is the origianl postion index + Hap length -1
  for (i in chr){
    PoS[[i]][,2]<-PoS[[i]][,2]+nHap-1
  }

  # Former
  # POS: current scanned snp positions (position, position+nHap-1)
  # POS.CHR: all snp position index in current chr (from full map)
  Extension<-function(Former, POS, POS.CHR, GEN, YFIX, KK, method, test, p.threshold, RR, PAR){
    WL<-list(NULL,NULL)
    #test.HAP resutls (wald or lrt) of scanned position
    FORMER<-Former[[test]]
    LATER<-FORMER
    # if not significant, return scanned position index range, scanned result, and NULL WL
    # if significant, continue extension
    if (FORMER[3]<p.threshold){
      iter<-0
      RE.WL<-list()
      # current chromsome,
      R.POS.CHR<-range(POS.CHR)
      # LATER/FORMER: statistics, left tail probability(log), p value, number of haplotype identified
      # loop end rule: LATER[1]<=FORMER[1] or LATER[2]>FORMER[2] or LATER[4]<=FORMER[4]
      # means: extension will end when this effect of haplotype is not increase anymore or numbers of haplotype identified less or equal than former
      while ((LATER[1]>FORMER[1] && LATER[2]<=FORMER[2] && LATER[4]>FORMER[4]) | iter==0){
        if (iter>0){
          FORMER<-LATER
          POS<-Pos
        }
        POS.T<-list()
        R.POS<-range(POS)
        # extent one bp on the left and right, not overflow the map range limit
        POS.ALL<-c(max(R.POS[1]-1,R.POS.CHR[1]):min(R.POS[2]+1,R.POS.CHR[2]))
        # left one bp and right one bp
        POS.LR<-list(c(min(POS.ALL):R.POS[2]),c(R.POS[1]:max(POS.ALL)))
        # the length of haplotype of left extension and right extension
        POS.LT<-lengths(POS.LR)
        # find the which direction is shorter
        # 1: left or equal
        # 2: right
        POS.MIN.LT<-which.min(POS.LT)
        # list(the shorter one, another one)
        POS.T<-list(POS.LR[[POS.MIN.LT]],POS.LR[[3-POS.MIN.LT]])

        later<-list()
        for (i in 1:2){
          # input numerical genotype matrix
          #GEN.x<-lapply(GK,function(x)x[,POS.T[[i]]])
          # position index range of haplotype
          POS.SE<-range(POS.T[[i]])

          POS.MAT<-which(RR[[test]][,1]==POS.SE[1] & RR[[test]][,2]==POS.SE[2])
          if (identical(POS.MAT,integer(0))){
            # test haplotype
            R.LATER<-test.HAP(t(GEN[POS.T[[i]],-(1:2)]),YFIX,KK,method,PAR)
            # add haplotype position index range to the test result
            T.LATER<-lapply(R.LATER[[1]],function(x)c(POS.SE,x))
          }else{
            T.LATER<-lapply(RR,function(x)x[POS.MAT[1],])
          }
          RE.WL[[2*iter+i]]<-T.LATER
          # test result
          later[[i]]<-T.LATER[[test]][-(1:2)]
        }
        # combine left and right extension result together
        Later<-do.call('rbind',later)
        # select most significant one according to left tail probability(log) and statistics, respectively.
        POS.PT<-c(which.min(Later[,2]),which.max(Later[,1]))
        # if statistic level is same, select the one with less snps (shorter haplotype)
        Pos.LR<-max(POS.PT[1],POS.PT[2-abs(diff(POS.LT))])
        Pos<-POS.T[[Pos.LR]]
        LATER<-Later[Pos.LR,]

        iter<-iter+1
      }
      # left direction extension result (all iteration)
      WL[[1]]<-do.call('rbind',lapply(RE.WL,function(x)x[[1]]))
      # right direction extension result (all iteration)
      WL[[2]]<-do.call('rbind',lapply(RE.WL,function(x)x[[2]]))
    }
    return(list(POS,FORMER,WL))
  }

  # if PAR is not provided, calculate it
  if (is.null(PAR)){
    PAR<-TEST.SCAN(YFIX,NULL,KIN,method,NULL)
  }

  n.test<-length(TEST)
  Final<-list(list(),list())
  RRS<-list()
  # scan and test haplotype position one by one, chr by chr
  for (i in chr){
    POS.CHR<-Pos.Chr[[i]]
    for (j in 1:nrow(PoS[[i]])){
      # scan position range
      pos<-PoS[[i]][j,1]:PoS[[i]][j,2]
      # extract snp code in pos range for all samples
      #GEN.x<-lapply(GK,function(x)x[,pos])
      RR0.S<-RR0[[TEST[1]]]
      RR0.E<-RR0[[TEST[n.test]]]
      pos.r<-range(pos)
      pos.temp<-which(RR0.S[,2]==pos.r[1] & RR0.S[,3]==pos.r[2] & RR0.E[,2]==pos.r[1] & RR0.E[,3]==pos.r[2])
      print(pos.temp)
      if (is.integer(pos.temp)){
        # identify haplotype?
        print(paste0("Current chr is ", i, " , Current initial position is ", pos))
        R.F<-test.HAP(t(GEN[pos,-(1:2)]),YFIX,KIN,method,PAR)
        #MAP0: chr, position index range, position itself in the range
        MAP0<-c(MAP[pos.r[1],1],pos.r,MAP[pos.r,2])
        rrs<-lapply(R.F[[1]],function(x)c(MAP0,x))
        RRS[[pos.r[1]]]<-rrs
      }else{
        RRS[[pos.r[1]]]<-lapply(RR0,function(x)x[pos.temp,])
      }
      # Former is the RF[[1]], which contain wald test and lrt test result
      # It is the initial value for extension
      Former<-lapply(RRS[[pos.r[1]]],function(x)x[-(1:5)])
      print(Former)
      RR<-list(NULL,NULL)
      for (ij in TEST){
        # Extension output list:
        # POS: the haplotype position index of the last significant test result
        # FORMER: the last significant test result
        # WL: all iteration result
        # Note: if initial is not significant, it will output initial result
        Re<-Extension(Former, pos, POS.CHR, GEN, YFIX, KIN, method, ij, p.threshold, RR, PAR)
        # significant haplotype map information
        pos.detail<-MAP[Re[[1]],1:2]
        # chr, position index, true position in map file, test result
        Final[[ij]][[pos.r[1]]]<-c(pos.detail[1,1],range(Re[[1]]),range(pos.detail[,2]),Re[[2]])
        RR<-Re[[3]]
      }

    }
  }


  RR0<-list()
  FINAL<-list()
  RR0[[1]]<-do.call('rbind',lapply(RRS,function(x)x[[1]]))
  RR0[[2]]<-do.call('rbind',lapply(RRS,function(x)x[[2]]))
  for (i in TEST){
    FINAL[[i]]<-do.call('rbind',Final[[i]])
  }
  return(list(FINAL,RR0))
}
