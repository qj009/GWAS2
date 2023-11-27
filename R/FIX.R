#' SNP statistical test function: FIXED
#' @description
#' This function allows you to do statistical test for selected SNP based on FIXED method. It is usually used inside the TEST.SCAN() function.
#' @details
#' Additional details...
#'

#' @param z: the input one snp genotype matrix for all samples, dim: 1*n, n is the sample counts.The rows represent samples. The columns represent SNPs. If z is NULL, then it will calculate PAR.
#' @param YFIX: Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.
#' @param KIN: Kinship matrix. It can be obtained from KIN() function.

#' @param fn0: Initial parameters for association test.

#' @returns: the test result out is a three-element list:
#' 1. wald test reuslt: wald test statistic, wald test left tail probability(log), wald test P value;
#' 2. liklihood ratio test(lrt) result: statistics, left tail probability(log), P value;
#' 3. parameters: beta, sigma2, lambda*sigma2, gamma, standand error, lrt statistics, lrt P value,wald test statistic,wald test P value.


#' @keywords cats
#' @export
#' @examples
#' cat_function()

#
FIX<-function(z,YFIX,KIN,fn0){
  options(digits=22)
  loglike<-function(theta){
    if (n.var>1){
      lambda<-as.list(exp(theta))
      Kin<-Map("*", KIN, lambda)
      KK<-Reduce("+", Kin)
      qq<-eigen(KK)
      yu<-t(qq[[2]])%*%y
      xu<-t(qq[[2]])%*%x
    }else{
      qq<-QQ
      qq[[1]]<-exp(theta)*qq[[1]]
    }

    logdt<-sum(log(qq[[1]]+1))
    h<-1/(qq[[1]]+1)
    yy<-sum(yu*h*yu)
    lx<-t(xu)%*%diag(h)
    yx<-lx%*%yu
    xx<-lx%*%xu
    CLxx<-tryCatch(chol(xx), error=function(e) NULL)
    if (is.null(CLxx)){
      ev<-eigen(xx)[[1]]
      INVxx<-ginv(xx)
    }else{
      ev<-diag(CLxx)^2
      INVxx<-chol2inv(CLxx)
    }
    DETxx<-prod(ev[ev>0])
    loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%INVxx%*%yx)-0.5*log(DETxx)
    if (!is.finite(loglike)) loglike<--1e+10
    return(-loglike)
  }

  fix<-function(theta){
    if (length(theta)>1){
      lambda<-as.list(exp(theta))
      Kin<-Map("*", KIN, lambda)
      KK<-Reduce("+", Kin)
      qq<-eigen(KK)
      yu<-t(qq[[2]])%*%y
      xu<-t(qq[[2]])%*%x
    }else{
      qq<-QQ
      qq[[1]]<-exp(theta)*qq[[1]]
    }

    h<-1/(qq[[1]]+1)
    yy<-sum(yu*h*yu)
    lx<-t(xu)%*%diag(h)
    yx<-lx%*%yu
    xx<-lx%*%xu

    CLxx<-tryCatch(chol(xx), error=function(e) NULL)
    if (is.null(CLxx)){INVxx<-ginv(xx)}else{INVxx<-chol2inv(CLxx)    }

    beta<-INVxx%*%yx
    sigma2<-(yy-t(yx)%*%INVxx%*%yx)/(n-q)
    sigma2<-drop(sigma2)
    stderr<-INVxx*sigma2
    return(list(beta,stderr,sigma2))
  }


  n<-nrow(YFIX)
  x<-as.matrix(YFIX[,-1])
  y<-as.matrix(YFIX[,1])
  n.var<-length(KIN)
  theta<-rep(0,n.var)
  Low<-rep(-10,n.var)
  Up<--Low
  if (n.var==1){
    QQ<-eigen(KIN[[1]])
    yu<-t(QQ[[2]])%*%y
  }
  if (is.null(z)){
    q<-ncol(x)
    if (n.var==1){xu<-t(QQ[[2]])%*%x}
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=Low,upper=Up)
    fn0<-parm$value
    return(fn0)
  }else{
    qx<-ncol(x)
    qz<-ncol(z)
    x<-cbind(x,z)
    q<-ncol(x)
    if (n.var==1){xu<-t(QQ[[2]])%*%x}
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=Low,upper=Up)
    Theta<-parm$par
    lambda<-exp(parm$par)
    conv<-parm$convergence
    hess<-parm$hessian
    parmfix<-fix(Theta)
    beta<-parmfix[[1]]
    stderr<-parmfix[[2]]
    sigma2<-parmfix[[3]]
    sigma2g<-lambda*sigma2
    fn1<-parm$value

    print(fn0)
    print(fn1)

    non.pos<-1:qx
    g<-beta[-non.pos]
    g.err<-stderr[-non.pos,-non.pos]
    b<-beta[non.pos]
    b.err<-stderr[non.pos]
    CLg.err<-tryCatch(chol(g.err), error=function(e) NULL)
    if (is.null(CLg.err)){
      INVg.err<-ginv(g.err)
    }else{
      INVg.err<-chol2inv(CLg.err)
    }
    wald<-t(g)%*%INVg.err%*%g
    p.w<-as.numeric(pchisq(wald,qz,log.p=TRUE))
    p.w<-abs(p.w)
    p.wald<-1-pchisq(wald,qz)
    WALD<-c(wald,p.w,p.wald)

    lrt<-2*(fn0-fn1)
    p.l<-as.numeric(pchisq(lrt,length(g),log.p=TRUE))
    p.l<-abs(p.l)
    p.lrt<-1-pchisq(lrt,length(g))
    LRT<-c(lrt,p.l,p.lrt)
    parm<-rep(0,9)
    parm[3]<-1
    RESULT<-list(WALD,LRT,parm)
    return(RESULT)
  }
}
