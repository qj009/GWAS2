# RANDOM method initial parameter calculation function: internal function
# @description
# This function allows you to calculate the initial parameter for RANDOM function. It is usually used inside the TEST.SCAN() function.
# @details
# Additional details...
#

# @param YFIX Phenotype input matrix. The first column is target phenotype data. The rest columns are FIXED traits user want to put into the model. If no FIXED traits, put 1 in the second column.
# @param KIN Kinship matrix. It can be obtained from KIN() function.


# @returns Initial parameter fn0 which can be used in RANDOM function.

# @import MASS
# @import  stats

MIXED<-function(YFIX,KIN){

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
    delta<-qq[[1]]
    logdt<-sum(log(delta+1))
    h<-1/(delta+1)
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


  x<-as.matrix(YFIX[,-1])
  y<-as.matrix(YFIX[,1])
  q<-ncol(x)
  n<-nrow(x)
  n.var<-length(KIN)
  if (n.var==1){
    QQ<-eigen(KIN[[1]])
    uu<-QQ[[2]]
    yu<-t(uu)%*%y
    xu<-t(uu)%*%x
  }
  theta<-rep(0,n.var)
  Parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=rep(-50,length(KIN)),upper=rep(10,length(KIN)))
  ##lambda<-exp(parm$par)
  return(Parm)
}
