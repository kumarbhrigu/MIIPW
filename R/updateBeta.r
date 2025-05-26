#' @title internal function for updating beta through Fisher Scoring
#' @param y response value for GEE model
#' @param x model matrix for the GEE model
#' @param vfun variance function for the GLM
#' @param mu mu vector for the GLM
#' @param w weight matrix
#' @param D derivation of the inverse link function
#' @param Ralpha correlation matrix
#' @param beta vector of beta value for GEE model
#' @export updateBeta
updateBeta<-function(y,x,vfun,mu,w,D,Ralpha,beta){
  hess<-t(D%*%x)%*%vfun%*%(Ralpha)%*%vfun%*%D%*%w%*%x
  u<-t(D%*%x)%*%vfun%*%(Ralpha)%*%vfun%*%w%*%as.matrix(y-mu)
  Newbeta<-beta+solve(as.matrix(hess),u)
  rslt<-list()
  rslt$Newbeta<-Newbeta
  rslt$hess<-hess
  rslt
}
