#' @title internal function for sandwich estimator
#' @details arguments are required for obtaining Sandwich Estimator for variance matrix of regression coefficient of GEE model
#' @param y response value for GEE model
#' @param x model matrix for the GEE model
#' @param vfun variance function for the GLM
#' @param mu mu vector for the GLM
#' @param w weight matrix
#' @param D derivation of the inverse link function
#' @param Ralpha correlation matrix
#' @param beta vector of beta value for GEE model
#' @param hessmat hessian matrix
#' @param blockdiag vector containing the dim of block matrix for block diagonal matrix
#'
#' @export updateSandW
updateSandW<-function(y,x,vfun,mu,w,D,Ralpha,beta,hessmat,blockdiag){
  blockmat<-function(blen,xmat){
    index<-blen
    xmat1<-xmat
    matlist<-list()
    for(i in 1:(length(blen)-1)){
      xmat2<-xmat1[1:blen[i],1:blen[i]]
      matlist[[i]]<-xmat2
      xmat1<-xmat1[-c(1:blen[i]),-c(1:blen[i])]
    }
    matlist[[length(blen)]]<-xmat1
    bmat<-bdiag(matlist)
    bmat

  }
  Blockdiag<-blockmat(blen=blockdiag,xmat=as.matrix(y-mu)%*%t(as.matrix(y-mu)))

  sigmahat<-(y-mu)%*%t(y-mu)
  usand<-t(x)%*%t(D)%*%vfun%*%Ralpha%*%vfun%*%w%*%sigmahat%*%w%*%vfun%*%Ralpha%*%vfun%*%D%*%x
  hessmat<-as.matrix(hessmat)
  upSandwich<-t(solve(hessmat,usand))
  #upSandwich<-upSandwich%*%solve(hessmat)
  upSandwich<-t(solve(t(hessmat),upSandwich))
  upSandwich
}
