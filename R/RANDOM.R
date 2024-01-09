#SNP statistical test function: RANDOM, internal functions
#description
#This function allows you to do statistical test for selected SNP based on RANDOM method. It is usually used inside the TEST.SCAN() function.
#details
#Additional details...
#'

#param z: the input one snp genotype matrix for all samples, dim: 1*n, n is the sample counts.The rows represent samples. The columns represent SNPs.
#param YFIX: Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.
#param KIN: Kinship matrix. It can be obtained from KIN() function.
#param theta: Initial parameters for association test.

#returns: the test result out is a three-element list:
#1. wald test reuslt: wald test statistic, wald test left tail probability(log), wald test P value;
#2. liklihood ratio test(lrt) result: statistics, left tail probability(log), P value;
#3. parameters: beta, sigma2, lambda*sigma2, gamma (two elements for heterozygote), standand error(four elements for heterozygote), wald test statistic,wald test P value,lrt statistics, lrt P value,.



RANDOM<-function(z,YFIX,KIN,Theta){
  Loglike<-function(theta){
    xi<-exp(theta)
    tt<-zz*xi+diag(r)
    CLtt<-tryCatch(chol(tt), error=function(e) NULL)
    if (is.null(CLtt)){
      ev<-eigen(tt)[[1]]
      INVtt<-ginv(tt)
    }else{
      ev<-diag(CLtt)^2
      INVtt<-chol2inv(CLtt)
    }
    DETtt<-prod(ev[ev>0])

    logdt2<-log(DETtt)
    zm<-xi*zx%*%INVtt
    my<-xi*INVtt%*%zy
    yHy<-yy-t(zy)%*%my
    yHx<-yx-zm%*%zy
    xHx<-xx-zm%*%t(zx)
    CLxHx<-tryCatch(chol(xHx), error=function(e) NULL)
    if (is.null(CLxHx)){
      ev<-eigen(xHx)[[1]]
      INVxHx<-ginv(xHx)
    }
    else{
      ev<-diag(CLxHx)^2
      INVxHx<-chol2inv(CLxHx)
    }
    DETxHx<-prod(ev[ev>0])
    loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%INVxHx%*%yHx)-0.5*log(DETxHx)
    if (!is.finite(loglike)) loglike<--1e+10
    return(-loglike)
  }

  fixed<-function(xi){
    tmp0<-zz*xi+diag(r)
    CLtmp0<-tryCatch(chol(tmp0), error=function(e) NULL)
    if (is.null(CLtmp0)){
      INVtmp0<-ginv(tmp0)
    }else{
      INVtmp0<-chol2inv(CLtmp0)
    }
    tmp<-xi*INVtmp0
    zxm<-zx%*%tmp
    zzm<-zz%*%tmp
    mzy<-tmp%*%zy
    yHy<-yy-t(zy)%*%mzy
    yHx<-yx-zx%*%mzy
    xHx<-xx-zxm%*%t(zx)
    zHx<-zx-zxm%*%zz
    zHy<-zy-zzm%*%zy
    zHz<-zz-zzm%*%zz
    CLxHx<-tryCatch(chol(xHx), error=function(e) NULL)
    if (is.null(CLxHx)){
      INVxHx<-ginv(xHx)
    }else{
      INVxHx<-chol2inv(CLxHx)
    }
    beta<-INVxHx%*%yHx
    sigma2<-(yHy-t(yHx)%*%INVxHx%*%yHx)/(n-s)
    gamma<-xi*zHy-xi*t(zHx)%*%INVxHx%*%yHx
    var<-abs((xi*diag(r)-xi*zHz*xi)*as.numeric(sigma2))
    stderr<-var
    result<-list(gamma,stderr,beta,sigma2)
    return(result)
  }

  x<-as.matrix(YFIX[,-1])
  y<-as.matrix(YFIX[,1])
  n<-nrow(y)
  s<-ncol(x)
  n.var<-length(KIN)

  lambda<-exp(Theta)
  if (n.var>1){
    lambda<-as.list(lambda)
    Kin<-Map("*", KIN, lambda)
    KK<-Reduce("+", Kin)
    qq<-eigen(KK)
  }else{
    qq<-eigen(KIN[[1]])
    qq[[1]]<-lambda*qq[[1]]
  }

  yu<-t(qq[[2]])%*%y
  xu<-t(qq[[2]])%*%x
  h<-1/(qq[[1]]+1)
  yy<-sum(yu*h*yu)
  lx<-t(xu)%*%diag(h)
  yx<-lx%*%yu
  xx<-lx%*%xu

  r<-ncol(z)
  zu<-t(qq[[2]])%*%z
  lz<-t(zu)%*%diag(h)
  zy<-lz%*%yu
  zz<-lz%*%zu
  zx<-t(lz%*%xu)

  Low<-rep(-10,1)
  Up<--Low

  theta<-rep(0,1)
  par<-optim(par=theta,fn=Loglike,hessian = TRUE,method="L-BFGS-B",lower=Low,upper=Up)
  xi<-exp(par$par)
  conv<-par$convergence
  fn1<-par$value
  hess<-par$hessian
  parmfix<-fixed(xi)

  gamma<-parmfix[[1]]
  stderr<-parmfix[[2]]
  beta<-parmfix[[3]]
  sigma2<-parmfix[[4]]
  lambda<-xi
  sigma2g<-lambda*sigma2
  ##fn0<-Parm$value
  fn0<-Loglike(-Inf)
  lrt<-2*(fn0-fn1)
  p.lrt<-1-pchisq(lrt,1)
  p.l<-abs(as.numeric(pchisq(lrt,1,log.p=TRUE)))
  CLstderr<-tryCatch(chol(stderr), error=function(e) NULL)
  if (is.null(CLstderr)){
    INVstderr<-ginv(stderr)
  }else{
    INVstderr<-chol2inv(CLstderr)
  }
  wald<-t(gamma)%*%INVstderr%*%gamma
  p.wald<-1-pchisq(wald,length(gamma))
  p.w<-abs(as.numeric(pchisq(wald,length(gamma),log.p=TRUE)))

  parm<-c(beta,sigma2,sigma2g,gamma,stderr,wald,p.wald,lrt,p.lrt)
  RESULT<-list(c(wald,p.w,p.wald),c(lrt,p.l,p.lrt),parm)
  return(RESULT)
}
