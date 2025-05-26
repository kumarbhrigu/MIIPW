#' QICmiipw
#' @title
#' Model Selection criteria QIC
#' @description It provides model selection criteria such as quasi-likelihood under the independence model criterion (QIC), an approximation to QIC under large sample i.e QICu and quasi likelihood
#' @param model.R fitted object obtained from GEE model \code{MeanScore,SIPW,AIPW,miSIPW,miAIPW} with correlation struture other than "independent"
#' @param model.indep same fitted object as in \code{model.indep} with "independent" correlation struture
#' @param family currently we have inlcuded
#'  "poisson","binomial","gaussian"
#' @return returns a list containing \code{QIC,QICu,Quasi likelihood}
#' @export
#' @references Pan, Wei. "Akaike's information criterion in generalized estimating equations." Biometrics 57.1 (2001): 120-125.
#' @examples
#' \dontrun{
#'  ##
#'  formula<-C6kine~ActivinRIB+ActivinRIIA+ActivinRIIAB+Adiponectin+AgRP+ALCAM
#'  pMat<-mice::make.predictorMatrix(srdata1[names(srdata1)%in%all.vars(formula)])
#'  m1<-MeanScore(data=srdata1,
#'              formula<-formula,id='ID',
#'              visit='Visit',family='gaussian',init.beta = NULL,
#'              init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
#'              corstr = 'exchangeable',maxit=50,m=2,pMat=pMat)
#'  m11<-MeanScore(data=srdata1,
#'              formula<-formula,id='ID',
#'              visit='Visit',family='gaussian',init.beta = NULL,
#'              init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
#'             corstr = 'independent',maxit=50,m=2,pMat=pMat)
#' QICmiipw(model.R=m1,model.indep=m11,family="gaussian")
#' ##
#' }
#'
QICmiipw <- function(model.R, model.indep,family) {
  fam<-family
  mu.R <- model.R$mu
  y <- model.R$y
  switch(family,
         poisson={quasiR <- sum((y*log(mu.R)) - mu.R)},
         gaussian={
           scale<-model.R$phi
           quasiR <- -1/2*sum((y-mu.R)^2/scale) },
         binomial={
           quasiR <- sum(y*log(mu.R/(1-mu.R))+log(1-mu.R))
         },stop("currently supported type poisson, gaussian, binomial ")
  )
  AIinv <- ginv(as.matrix(model.indep$betaSand))# Omega-hat(I) via Moore-Penrose
    vr <- model.R$betaSand
  traceR <- sum(diag(AIinv %*% vr))
  px <- length(mu.R)

  # QIC
  QIC <- (-2)*quasiR + 2*traceR
  QICu <- (-2)*quasiR + 2*px
  result <- c(QIC, QICu, quasiR)
  names(result) <- c('QIC', 'QICu', 'Quasi Lik')
  result}
