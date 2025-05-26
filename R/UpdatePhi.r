#' @title internal function for updating scale parameter
#' @param y response value for GEE model
#' @param x model matrix for the GEE model
#' @param vfun variance function for the GLM
#' @param mu mu vector for the GLM
#' @param w weight matrix
#' @export UpdatePhi
UpdatePhi<-function(y,x,vfun,mu,w){
  ymu<-(as.matrix(y-mu,ncol=1,nrow=length(y),byrow=T))
  res<-vfun%*%(w%*%(ymu))
  Newphi<-sum(res*res)/(length(y)-ncol(x))
  Newphi
}
