#' lmemeanscore
#' @title
#' Fits a marginal model using meanscore
#' @description provides meanscore estimates of parameters for semiparametric marginal model
#' of response variable of interest. The augmented terms are estimated by using multiple imputation model.
#' @details It uses the mean score method to reduce the bias
#' due to missing values in response model for longitudinal data.The response variable \eqn{\mathbf{Y}} is related to the coariates as \eqn{g(\mu)=\mathbf{X}\beta}, where \code{g} is the link function for the glm. The estimating equation is
#' \deqn{\sum_{i=1}^{n}\sum_{j=t_1}^{t_k}\int_{b_i}(\delta_{ij}S(Y_{ij},\mathbf{X}_{ij})+(1-\delta_{ij})\phi(\mathbf{V}_{ij},b_i;\psi))db_i=0}
#' where \eqn{\delta_{ij}=1} if there is missing value in response and 0 otherwise,
#' \eqn{\mathbf{X}} is fully observed all subjects and
#' \eqn{\mathbf{V}_{ij}=(\mathbf{X}_{ij},A_{ij})}. The missing score function values due to incomplete data are estimated
#' using an imputation model through mice which we have considered as \deqn{Y_{ij}|b_i\sim N(\mathbf{V}_{ij}\gamma+b_i,\sigma)\; ; b_i\sim N(0,\sigma_{miss})}
#' through multiple imputation.
#' @param data longitudinal data with each subject specified discretely
#' @param M number of imputation to be used in the estimation of augmentation term
#' @param id cloumn names which shows identification number for each subject
#' @param analysis.model A formula to be used as analysis model
#' @param imp.model For for missing response imputation, which consider subject specific random intercept
#' @param qpoints Number of quadrature points to be used while evaluating the numerical integration
#' @param psiCov working model parameter
#' @param psi working model parameter
#' @param sigma working model parameter
#' @param sigmaMiss working model parameter
#' @param dist distribution for imputation model. Currently available options are Gaussian, Binomial
#' @param link Link function for the mean
#' @param conv convergence tolerance
#' @param maxiter  maximum number of iteration
#' @param maxpiinv maximum value pi can take
#' @param se Logical for Asymptotic SE for regression coefficient of the regression model.
#' @param verbose logical argument
#' @return A list of objects containing the following objects
#' \describe{
#'   \item{Call}{details about arguments passed in the function}
#'   \item{nr.conv}{logical for checking convergence in Newton Raphson algorithm}
#'   \item{nr.iter}{number of iteration required}
#'   \item{nr.diff}{absolute difference for roots of Newton Raphson algorithm}
#'   \item{beta}{estimated regression coefficient for the analysis model}
#'   \item{var.beta}{Asymptotic SE for beta}
#' }
#' @import lme4
#' @import spatstat
#' @export
#' @references Wang, C. Y., Shen-Ming Lee, and Edward C. Chao. "Numerical equivalence of imputing scores and weighted estimators in regression analysis with missing covariates." Biostatistics 8.2 (2007): 468-473.
#' @references Seaman, Shaun R., and Stijn Vansteelandt. "Introduction to double robust methods for incomplete data." Statistical science: a review journal of the Institute of Mathematical Statistics 33.2 (2018): 184.
#' @references Vansteelandt, Stijn, James Carpenter, and Michael G. Kenward. "Analysis of incomplete data using inverse probability weighting and doubly robust estimators." Methodology: European Journal of Research Methods for the Behavioral and Social Sciences 6.1 (2010): 37.
#' @examples
#'  \dontrun{
#' ##
#' library(JMbayes2)
#' library(lme4)
#' library(insight)
#' library(numDeriv)
#' library(stats)
#' lmer(log(alkaline)~drug+age+year+(1|id),data=na.omit(pbc2))
#' data1<-pbc2
#' data1$alkaline<-log(data1$alkaline)
#' names(pbc2)
#' apply(pbc2,2,function(x){sum(is.na(x))})
#' r.ij<-ifelse(is.na(data1$alkaline)==T,0,1)
#' data1<-cbind.data.frame(data1,r.ij)
#' data1$drug<-factor(data1$drug,levels=c("placebo","D-penicil"),labels = c(0,1))
#' data1$sex<-factor(data1$sex,levels=c('male','female'),labels=c(1,0))
#' data1$drug<-as.numeric(as.character(data1$drug))
#' data1$sex<-as.numeric(as.character(data1$sex))
#' model.y<-lmer(alkaline~year+age+sex+drug+serBilir+(1|id),data=na.omit(data1))
#' psi<-model.y@beta
#' sigma<-get_variance_residual(model.y)
#' sigmaMiss<-get_variance(model.y)$var.random
#' m11<-lmemeanscore(data=data1,id='id',
#' analysis.model = alkaline~year,
#' imp.model = ~year+age+sex+drug+serBilir+(1|id),
#' psiCov = vcov(model.y),psi=psi,
#' sigma=sigma,sigmaMiss=sigmaMiss,dist='gaussian',link='identity',qpoints = 4,
#' maxiter = 200)
#' m11
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso   \link[MIIPW]{SIPW},\link[MIIPW]{miSIPW},\link[MIIPW]{miAIPW}
lmemeanscore<-function(data,M=5,id,analysis.model,imp.model,
                       qpoints=4,psiCov,psi,sigma=NULL,
                       sigmaMiss,dist,link,
                       conv=.0001,maxiter,
                       maxpiinv=-1,se=TRUE,verbose=FALSE){
  #library(spatstat)
  arg_checks <- as.list(match.call())[-1]
  arg_checks$analysis.model<-deparse(analysis.model)
  #arg_checks$wgt.model<-deparse(wgt.model)
  arg_checks$imp.model<-deparse(imp.model)
  if(is.null(data)){
    stop("data cannot be null, provide dataset which is used
         in the response model as well as in the working model")
  }
  if(is.null(id)){
    stop("id cannot be null, provide a longitudinal data
         with individual identificatio number")
  }
  if(is.null(analysis.model)){
    stop("provide an analysis model")
  }

  if(is.null(imp.model)){
    stop("provide an imputation model")
  }

  if(is.null(psi)||is.null(psiCov)||is.null(sigmaMiss)){
    stop("provide working model parameters")
  }
  if(is.null(dist)||is.null(link)){
    stop("provide distribution and link function compatible with the data provided")
  }
  if(id%in%names(data)){
    names(data)[which(names(data)==id)]<-'id'
  }
  #match(c('formula', 'data', 'weights'), names(mf),0L)
  data<-data
  dist<-dist
  link<-link
  gammaCov<-psiCov
  gamma<-psi
  sigma<-sigma
  phi<-sigmaMiss

  analysis.var<-all.vars(analysis.model)
  imp.var<-all.vars(nobars(imp.model))

  if(sum(analysis.var%in%names(data))!=length(analysis.var)){
    stop('variables included in the models used should be in the dataset')
  }



  if(sum(imp.var%in%names(data))!=length(imp.var)){
    stop('variables included in the models used should be in the dataset')
  }

  y<-data[,analysis.var[1]]
  x<-as.matrix(cbind(1,data[,c(analysis.var[-1])]))
  z.y<-data[,c(imp.var)]
  n<-nrow(data)
  m<-length(unique(data$id))
  nj<-as.data.frame(table(data$id))$Freq
  # ASSIGN VARIABLES
  id<-rep(1:m,times=nj)
  r<-1*!is.na(y)
  covars.y<-cbind(intercept=rep(1,n),z.y);covars.y<-as.matrix(covars.y)
  p<-c(ncol(x),0,ncol(covars.y))
  nj.obs<-aggregate(r~id,FUN=sum,data=data.frame(r,id))$r
  n.obs<-sum(r)
  sigma<-ifelse(sigma>=.0001,sigma,.0001);  phi<-ifelse(phi>=.0001,phi,.0001)
  gammalist<-mvrnorm(n=M,mu=gamma,Sigma=gammaCov)
  # ESTIMATE INTEGRALS FROM ESTIMATING EQUATIONS USING GAUSS-HERMITE QUADRATURE
  if(dist=='gaussian'){
    obs.resid<-ifelse(r,(y-covars.y%*%gamma),0)
    obs.sq.resid<-obs.resid^2
    int.b<-(nj.obs/sigma+1/phi)^(-1)*aggregate(obs.resid~id,FUN=sum,data=data.frame(obs.resid,id))$obs.resid/sigma
    int.nulist<-list()
    for(i in 1:M){
      int.nulist[[i]]<-c(covars.y%*%as.vector(gammalist[i,]))+int.b[id]
    }
    int.nu<-Reduce('+',int.nulist)/M
  }

  if(verbose){print('Posterior Expectation of Nu'); print(summary(int.nu))}
  trim=NULL

  # NEWTON-RAPHSON TO ESTIMATE BETA FROM ESTIMATING EQUATIONS
  beta<-glm(y ~.-1, family=dist, data=cbind.data.frame(y,x))$coefficients
  nr.iter<-0;
  diff.nr<-100;
  while((max(abs(diff.nr))>conv)&(nr.iter<maxiter)){

    if(link=='identity'){
      mu.beta<-c(x%*%beta)
      dmu.beta<-x
      d2mu.beta<-matrix(0,n,p)
    }

    gij<-dmu.beta*c(ifelse(r,y-int.nu,0))+dmu.beta*int.nu-dmu.beta*mu.beta
    sm.beta<-colSums(gij)
    dsm.beta<-t(x)%*%(d2mu.beta*c(ifelse(r,y-int.nu,0))+int.nu-mu.beta)-t(dmu.beta)%*%dmu.beta
    beta.new<-as.vector(beta-solve(dsm.beta)%*%sm.beta)
    diff.nr<-beta.new-beta
    beta=beta.new; remove(beta.new); nr.iter<-nr.iter+1;
    if(verbose){print('Newton-Raphson'); print(c(nr.iter,beta))}
  }
  if(max(abs(diff.nr))>conv){print(paste0('WARNING: Newton-Raphson algorithm did not converge. Maximum absolute change in parameter estimates at the last iteration was ',max(abs(diff.nr))))}

  SE<-function(beta,gamma,sigma,phi){
    g<-function(gamma,sigma,phi){
      sigma<-ifelse(sigma>=.0001,sigma,.0001);  phi<-ifelse(phi>=.0001,phi,.0001)
      if(dist=='gaussian'){	# link: identity
        obs.resid<-ifelse(r,(y-covars.y%*%gamma),0)
        obs.sq.resid<-obs.resid^2
        int.b<-(nj.obs/sigma+1/phi)^(-1)*aggregate(obs.resid~id,FUN=sum,data=data.frame(obs.resid,id))$obs.resid/sigma
        int.nulist<-list()
        for(i in 1:M){
          int.nulist[[i]]<-c(covars.y%*%as.vector(gammalist[i,]))+int.b[id]
        }
        int.nu<-Reduce('+',int.nulist)/M
      }
      if(link=='identity'){
        mu.beta<-c(x%*%beta)
        dmu.beta<-x
        d2mu.beta<-matrix(0,n,p)
      }
      if(link=='logit'){
        mu.beta<-1/(exp(-c(x%*%beta))+1)
        dmu.beta<-x*exp(-c(x%*%beta))/(exp(-c(x%*%beta))+1)^2
        d2mu.beta<-x*exp(-c(x%*%beta))/(exp(-c(x%*%beta))+1)^3
      }
      gij<-dmu.beta*c(ifelse(r,y-int.nu,0))+dmu.beta*int.nu-dmu.beta*mu.beta
      return(apply(gij,2,sum)/m)
    }
    env.g=new.env()
    assign("gamma", gamma, envir = env.g)
    assign("sigma", sigma, envir = env.g)
    assign("phi", phi, envir = env.g)
    if(dist=='gaussian'){
      d1.gij<-attr(numericDeriv(quote(g(gamma,sigma,phi)), c('gamma','sigma','phi'), env.g),'gradient')     # 1st derivative of E[within-cluster-sum(g_ij)]
      colnames(d1.gij)<-c(paste0(rep('gamma',p[3]),1:p[3]),'sigma','phi')
    }

    if(verbose){print('Derivative of gij')}

    # derivatives for l_j(gamma,phi)
    if(dist=='gaussian'){
      l.gamma.phi.sigma.each<-function(gamma,phi,sigma){
        phi<-ifelse(phi>=.0001,phi,.0001)
        y.pdf<-function(b){    # binary/Bernoulli distribution (logit link)
          y.pdf.each<-ifelse(r,exp(-(y-c(covars.y%*%gamma)-b)^2/(2*sigma))*(1/(sqrt(2*pi*sigma))),1)
          return(aggregate(y.pdf.each~id,FUN=prod,data=data.frame(id,y.pdf.each))$y.pdf.each)
        }
        return(log(gauss.hermite(y.pdf,mu=0,sd=sqrt(phi),order=qpoints)))
      }
      l.gamma.phi.sigma<-function(parm){
        gamma=parm[1:(p[3])]; phi=parm[p[3]+1]; sigma=parm[p[3]+2]
        return(mean(l.gamma.phi.sigma.each(gamma,phi,sigma)))
      }
      env.l.gamma.phi.sigma=new.env()
      assign("gamma", gamma, envir = env.l.gamma.phi.sigma)
      assign("phi", phi, envir = env.l.gamma.phi.sigma)
      assign("sigma", sigma, envir = env.l.gamma.phi.sigma)
      d1l.gamma.phi<-attr(numericDeriv(quote(l.gamma.phi.sigma.each(gamma,phi,sigma)), c('gamma','phi','sigma'), env.l.gamma.phi.sigma),'gradient')

      d2l.gamma.phi<-hessian(l.gamma.phi.sigma,x=c(gamma,phi,sigma))
      colnames(d1l.gamma.phi)<-c(paste0(rep('gamma',p[3]),1:p[3]),'phi','sigma'); colnames(d2l.gamma.phi)=colnames(d1l.gamma.phi); rownames(d2l.gamma.phi)=colnames(d1l.gamma.phi)
    }
    if(verbose){print('1st Derivative of l(gamma,phi)')}
    if(verbose){print('2nd Derivative of l(gamma,phi)')}
    gj<-aggregate(gij~id,FUN=sum,data=data.frame(id,gij))[,2:(p[1]+1)]
    inside<-solve(dsm.beta)%*%(t(gj)-d1.gij[,grepl('gamma',colnames(d1.gij))|colnames(d1.gij)=='sigma'|colnames(d1.gij)=='phi']%*%solve(d2l.gamma.phi)%*%t(d1l.gamma.phi))
    var.beta<-inside%*%t(inside)
    colnames(var.beta)<-c(paste0(rep('beta',p[1]),0:(p[1]-1))); rownames(var.beta)=colnames(var.beta)
    return(var.beta)
  }
  if(se==TRUE){var.beta<-SE(beta,gamma,sigma,phi)}
  if(se==FALSE){var.beta=NULL}
  fitvalues<-c(as.matrix(x)%*%as.matrix(beta,ncol=1))
  gammalist<-gammalist
  #result=list(Call=arg_checks,nr.conv=(max(abs(diff.nr))<=conv),
  #            nr.iter=nr.iter,nr.diff=max(abs(diff.nr)),beta=beta,var.beta=var.beta,fitvalues=fitvalues,gammalist=gammalist)


  result=list(Method='lmemeanscore',Call=arg_checks,nr.conv=(max(abs(diff.nr))<=conv),nr.iter=nr.iter,nr.diff=max(abs(diff.nr)),beta=beta,var.beta=var.beta,fitvalues=fitvalues,gammalist=gammalist)

  return(result)

}

utils::globalVariables(c('data1','aggregate','gauss.hermite','numericDeriv','hessian'))
