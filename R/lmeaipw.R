
#' lmeaipw
#' @title
#' Fits a marginal model using AIPW
#' @description provides augmented inverse probability weighted estimates of parameters for semiparametric marginal model
#' of response variable of interest. The augmented terms are estimated by using multiple imputation model.
#' @details It uses the augmented inverse probability weighted method to reduce the bias
#' due to missing values in response model for longitudinal data. The response variable \eqn{\mathbf{Y}} is related to the coariates as \eqn{g(\mu)=\mathbf{X}\beta}, where \code{g} is the link function for the glm. The estimating equation is
#' \deqn{\sum_{i=1}^{n}\sum_{j=t_1}^{t_k}\int_{a_i}\int_{b_i}(\frac{\delta_{ij}}{\hat\pi_{ij}(a_i)}S(Y_{ij},\mathbf{X}_{ij};\beta)+(1-\frac{\delta_{ij}}{\hat\pi_{ij}(a_i)})\phi(\mathbf{V}_{ij},b_i;\psi))da_idb_i=0}
#' where \eqn{\delta_{ij}=1} if there is missing value in the response and 0 otherwise,
#' \eqn{\mathbf{X}} is fully observed all subjects,
#' where \eqn{\mathbf{V}_{ij}=(\mathbf{X}_{ij},A_{ij})}. The missing score function values due to incomplete data are estimated
#' using an imputation model through FCS (here we have considered a mixed effect model) which we have considered as \eqn{\phi(\mathbf{V}_{ij}=\mathbf{v}_{ij}))}. The estimated value \eqn{\hat{\phi}(\mathbf{V}_{ij}=\mathbf{v}_{ij}))} is obtained
#' through multiple imputation. The working model for imputation of missing response is \deqn{Y_{ij}|b_i\sim N(\mathbf{V}_{ij}\psi+b_i,\sigma)\; ; b_i\sim N(0,\sigma_{miss})} and
#' for the missing data probability \deqn{Logit(P(\delta_{ij}=1|\mathbf{V}_{ij}\nu+a_i))\;;a_i\sim N(0,\sigma_R)}
#' @param data longitudinal data with each subject specified discretely
#' @param M number of imputation to be used in the estimation of augmentation term
#' @param id cloumn names which shows identification number for each subject
#' @param analysis.model A formula to be used as analysis model
#' @param wgt.model Formula for weight model, which consider subject specific random intercept
#' @param imp.model For for missing response imputation, which consider subject specific random intercept
#' @param qpoints Number of quadrature points to be used while evaluating the numerical integration
#' @param psiCov working model parameter
#' @param nu working model parameter
#' @param psi working model parameter
#' @param sigma working model parameter
#' @param sigmaMiss working model parameter
#' @param sigmaR working model parameter
#' @param dist distribution for imputation model. Currently available options are Gaussian, Binomial
#' @param link Link function for the mean
#' @param conv convergence tolerance
#' @param maxiter maximum number of iteration
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
#' r.ij~year+age+sex+drug+serBilir+(1|id)
#' model.r<-glmer(r.ij~year+age+sex+drug+serBilir+(1|id),family=binomial(link='logit'),data=data1)
#' model.y<-lmer(alkaline~year+age+sex+drug+serBilir+(1|id),data=na.omit(data1))
#' nu<-model.r@beta
#' psi<-model.y@beta
#' sigma<-get_variance_residual(model.y)
#' sigmaR<-get_variance(model.r)$var.random
#' sigmaMiss<-get_variance(model.y)$var.random
#' m11<-lmeaipw(data=data1,id='id',
#' analysis.model = alkaline~year,
#' wgt.model=~year+age+sex+drug+serBilir+(1|id),
#' imp.model = ~year+age+sex+drug+serBilir+(1|id),
#' psiCov = vcov(model.y),nu=nu,psi=psi,
#' sigma=sigma,sigmaMiss=sigmaMiss,sigmaR=sigmaR,dist='gaussian',link='identity',
#' maxiter = 200)
#' m11
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso   \link[MIIPW]{SIPW},\link[MIIPW]{miSIPW},\link[MIIPW]{miAIPW}
lmeaipw<-function(data,M=5,id,
                  analysis.model,wgt.model,imp.model,
                  qpoints=4,
                  psiCov,nu,psi,
                  sigma=NULL,sigmaMiss,sigmaR,dist,link,
                  conv=.0001,maxiter,maxpiinv=-1,
                  se=TRUE,verbose=FALSE){
  # assumes independent normally-distributed random intercepts only, two-level data
  #library(spatstat)

  arg_checks <- as.list(match.call())[-1]
  arg_checks$analysis.model <-deparse(analysis.model)
  arg_checks$wgt.model<-deparse(wgt.model)
  arg_checks$imp.model<-deparse(imp.model)
  #arg_checks$psiCov<-

  if(is.null(data)){
    stop('data cannot be null, provide dataset which is used
         in the response model as well as in the working model')
  }
  if(is.null(id)){
    stop('id cannot be null, provide a longitudinal data
         with individual identification number')
  }
  if(is.null(analysis.model)){
    stop('provide an analysis model')
  }

  if(is.null(imp.model)){
    stop('provide an imputation model')
  }
  if(is.null(wgt.model)){
    stop('provide a weight model')
  }
  if(is.null(psiCov)||is.null(nu)||is.null(psi)||is.null(sigmaR)||is.null(sigmaMiss)){
    stop('provide working model parameters')
  }
  if(is.null(dist)||is.null(link)){
    stop('provide distribution and link function compatible with the data provided')
  }
  if(id%in%names(data)){
    names(data)[which(names(data)==id)]<-'id'
  }
  #match(c('formula', 'data', 'weights'), names(mf),0L)
  data<-data
  dist<-dist
  link<-link
  gammaCov<-psiCov
  alpha<-nu
  gamma<-psi
  sigma<-sigma
  tau<-sigmaR
  phi<-sigmaMiss

  analysis.var<-all.vars(analysis.model)
  wgt.var<-all.vars(nobars(wgt.model))
  imp.var<-all.vars(nobars(imp.model))

  if(sum(analysis.var%in%names(data))!=length(analysis.var)){
    stop('variables included in the models used should be in the dataset')
  }

  if(sum(wgt.var%in%names(data))!=length(wgt.var)){
    stop('variables included in the models used should be in the dataset')
  }

  if(sum(imp.var%in%names(data))!=length(imp.var)){
    stop('variables included in the models used should be in the dataset')
  }


  y<-data[,analysis.var[1]]
  x<-as.matrix(cbind(1,data[,c(analysis.var[-1])]))
  z.r<-data[,c(wgt.var)]
  z.y<-data[,c(imp.var)]
  n<-nrow(data)
  m<-length(unique(data$id))
  nj<-as.data.frame(table(data$id))$Freq
  # ASSIGN VARIABLES
  id<-rep(1:m,times=nj)
  r<-1*!is.na(y)
  covars.r<-cbind(intercept=rep(1,n),z.r);covars.r<-as.matrix(covars.r)
  covars.y<-cbind(intercept=rep(1,n),z.y);covars.y<-as.matrix(covars.y)
  p<-c(ncol(x),ncol(covars.r),ncol(covars.y))
  nj.obs<-aggregate(r~id,FUN=sum,data=data.frame(r,id))$r
  n.obs<-sum(r)
  sigma<-ifelse(sigma>=.0001,sigma,.0001); tau<-ifelse(tau>=.0001,tau,.0001); phi<-ifelse(phi>=.0001,phi,.0001)
  gammalist<-mvrnorm(n=M,mu=gamma,Sigma=gammaCov)
  if(dist=='gaussian'){	# link: identity
    obs.resid<-ifelse(r,(y-covars.y%*%gamma),0)
    obs.sq.resid<-obs.resid^2
    int.b<-(nj.obs/sigma+1/phi)^(-1)*aggregate(obs.resid~id,FUN=sum,data=data.frame(obs.resid,id))$obs.resid/sigma  # posterior expectation for b
    int.nulist<-list()
    for(i in 1:M){
      int.nulist[[i]]<-c(covars.y%*%as.vector(gammalist[i,]))+int.b[id]
    }
    int.nu<-Reduce('+',int.nulist)/M
  }
  if(verbose){print('Posterior Expectation of Nu'); print(summary(int.nu))}
  r.pdf<-function(a){
    r.pdf.each<-exp(r*(c(covars.r%*%alpha)+a))/(1+exp(c(covars.r%*%alpha)+a))
    return(aggregate(r.pdf.each~id,FUN=prod,data=data.frame(id,r.pdf.each))$r.pdf.each)
  }
  r.pdf.overpi<-function(a){
    invpi<-exp(-c(covars.r%*%alpha)-a)+1
    return(invpi*r.pdf(a)[id])
  }
  int.piinv<-gauss.hermite(r.pdf.overpi,mu=0,sd=sqrt(tau),order=qpoints)/gauss.hermite(r.pdf,mu=0,sd=sqrt(tau),order=qpoints)[id]   # posterior expectation for 1/pi
  trim=NULL; if(maxpiinv>0){trim=sum(int.piinv>maxpiinv); int.piinv[int.piinv>maxpiinv]=maxpiinv}		# trim inverse of pi to avoid extremely influential values
  if(verbose){print('Posterior Expectation of Inverse of Pi'); print(summary(int.piinv))}
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
    gij<-dmu.beta*c(ifelse(r,y*int.piinv-int.nu*int.piinv,0))+dmu.beta*int.nu-dmu.beta*mu.beta
    sm.beta<-colSums(gij)
    dsm.beta<-t(x)%*%(d2mu.beta*c(ifelse(r,y*int.piinv-int.nu*int.piinv,0))+int.nu-mu.beta)-t(dmu.beta)%*%dmu.beta
    beta.new<-as.vector(beta-solve(dsm.beta)%*%sm.beta)
    diff.nr<-beta.new-beta
    beta<-beta.new; remove(beta.new); nr.iter<-nr.iter+1;
    if(verbose){print('Newton-Raphson'); print(c(nr.iter,beta))}
  }
  if(max(abs(diff.nr))>conv){print(paste0('WARNING: Newton-Raphson algorithm did not converge. Maximum absolute change in parameter estimates at the last iteration was ',max(abs(diff.nr))))}
  SE<-function(beta,alpha,gamma,sigma,tau,phi){
    #library(spatstat)
    # FUNCTION TO ESTIMATE ASYMPTOTIC COVARIANCE MATRIX
    # derivatives for g_ij
    g<-function(alpha,gamma,sigma,tau,phi){

      sigma<-ifelse(sigma>=.0001,sigma,.0001); tau=ifelse(tau>=.0001,tau,.0001); phi=ifelse(phi>=.0001,phi,.0001)

      # ESTIMATE INTEGRALS FROM ESTIMATING EQUATIONS USING GAUSS-HERMITE QUADRATURE
      if(dist=='gaussian'){	# link: identity
        obs.resid<-ifelse(r,(y-covars.y%*%gamma),0)
        obs.sq.resid<-obs.resid^2
        int.b<-(nj.obs/sigma+1/phi)^(-1)*aggregate(obs.resid~id,FUN=sum,data=data.frame(obs.resid,id))$obs.resid/sigma  # posterior expectation for b
        int.nulist<-list()
        for(i in 1:M){
          int.nulist[[i]]<-c(covars.y%*%as.vector(gammalist[i,]))+int.b[id]
        }
        int.nu=Reduce('+',int.nulist)/M
      }

      if(dist=='binomial'){
        y.pdf<-function(b){
          y.pdf.each<-ifelse(r,exp(y*(c(covars.y%*%gamma)+b))/(1+exp(c(covars.y%*%gamma)+b)),1)
          return(aggregate(y.pdf.each~id,FUN=prod,data=data.frame(id,y.pdf.each))$y.pdf.each)
        }
        y.pdf.nu<-function(b){
          nu<-1/(exp(-c(covars.y%*%gamma)-b)+1)
          return(nu*y.pdf(b)[id])
        }
        int.nu<-gauss.hermite(y.pdf.nu,mu=0,sd=sqrt(phi),order=qpoints)/gauss.hermite(y.pdf,mu=0,sd=sqrt(phi),order=qpoints)[id]
      }

      r.pdf<-function(a){
        r.pdf.each<-exp(r*(c(covars.r%*%alpha)+a))/(1+exp(c(covars.r%*%alpha)+a))
        return(aggregate(r.pdf.each~id,FUN=prod,data=data.frame(id,r.pdf.each))$r.pdf.each)
      }
      r.pdf.overpi<-function(a){
        invpi<-exp(-c(covars.r%*%alpha)-a)+1
        return(invpi*r.pdf(a)[id])
      }
      int.piinv<-gauss.hermite(r.pdf.overpi,mu=0,sd=sqrt(tau),order=qpoints)/gauss.hermite(r.pdf,mu=0,sd=sqrt(tau),order=qpoints)[id]


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

      gij<-dmu.beta*c(ifelse(r,y*int.piinv-int.nu*int.piinv,0))+dmu.beta*int.nu-dmu.beta*mu.beta
      return(apply(gij,2,sum)/m)
    }
    env.g<-new.env()
    assign("alpha", alpha, envir = env.g)
    assign("gamma", gamma, envir = env.g)
    assign("sigma", sigma, envir = env.g)
    assign("tau", tau, envir = env.g)
    assign("phi", phi, envir = env.g)
    if(dist=='gaussian'){
      d1.gij<-attr(numericDeriv(quote(g(alpha,gamma,sigma,tau,phi)), c('alpha','gamma','sigma','tau','phi'), env.g),'gradient')
      colnames(d1.gij)<-c(paste0(rep('alpha',p[2]),1:p[2]),paste0(rep('gamma',p[3]),1:p[3]),'sigma','tau','phi')
    }
    if(dist=='binomial'){
      d1.gij<-attr(numericDeriv(quote(g(alpha,gamma,sigma,tau,phi)), c('alpha','gamma','tau','phi'), env.g),'gradient')
      colnames(d1.gij)<-c(paste0(rep('alpha',p[2]),1:p[2]),paste0(rep('gamma',p[3]),1:p[3]),'tau','phi')
    }
    if(verbose){print('Derivative of gij')}

    l.alpha.tau.each<-function(alpha,tau){
      tau<-ifelse(tau>=.0001,tau,.0001)
      r.pdf<-function(a){
        r.pdf.each<-exp(r*(c(covars.r%*%alpha)+a))/(1+exp(c(covars.r%*%alpha)+a))
        return(aggregate(r.pdf.each~id,FUN=prod,data=data.frame(id,r.pdf.each))$r.pdf.each)
      }
      return(log(gauss.hermite(r.pdf,mu=0,sd=sqrt(tau),order=qpoints)))
    }
    l.alpha.tau<-function(parm){
      alpha<-parm[1:(p[2])]; tau=parm[p[2]+1]
      return(mean(l.alpha.tau.each(alpha,tau)))
    }
    env.l.alpha.tau<-new.env()
    assign("alpha", alpha, envir = env.l.alpha.tau)
    assign("tau", tau, envir = env.l.alpha.tau)
    d1l.alpha.tau<-attr(numericDeriv(quote(l.alpha.tau.each(alpha,tau)), c('alpha','tau'), env.l.alpha.tau),'gradient')   # 1st derivative of l_j(alpha,tau) for each cluster j=1,...,m
    if(verbose){print('1st Derivative of l(alpha,tau)')}

    d2l.alpha.tau<-hessian(l.alpha.tau,x=c(alpha,tau))
    colnames(d1l.alpha.tau)<-c(paste0(rep('alpha',p[2]),1:p[2]),'tau'); colnames(d2l.alpha.tau)=colnames(d1l.alpha.tau); rownames(d2l.alpha.tau)=colnames(d1l.alpha.tau)
    if(verbose){print('2nd Derivative of l(alpha,tau)')}
    if(dist=='gaussian'){
      l.gamma.phi.sigma.each<-function(gamma,phi,sigma){
        phi<-ifelse(phi>=.0001,phi,.0001)
        y.pdf<-function(b){
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
      d1l.gamma.phi=attr(numericDeriv(quote(l.gamma.phi.sigma.each(gamma,phi,sigma)), c('gamma','phi','sigma'), env.l.gamma.phi.sigma),'gradient')
      d2l.gamma.phi=hessian(l.gamma.phi.sigma,x=c(gamma,phi,sigma))
      colnames(d1l.gamma.phi)=c(paste0(rep('gamma',p[3]),1:p[3]),'phi','sigma'); colnames(d2l.gamma.phi)=colnames(d1l.gamma.phi); rownames(d2l.gamma.phi)=colnames(d1l.gamma.phi)
    }
    if(dist=='binomial'){
      l.gamma.phi.each<-function(gamma,phi){
        phi<-ifelse(phi>=.0001,phi,.0001)
        y.pdf<-function(b){
          y.pdf.each<-ifelse(r,exp(y*(c(covars.y%*%gamma)+b))/(1+exp(c(covars.y%*%gamma)+b)),1)
          return(aggregate(y.pdf.each~id,FUN=prod,data=data.frame(id,y.pdf.each))$y.pdf.each)
        }
        return(log(gauss.hermite(y.pdf,mu=0,sd=sqrt(phi),order=qpoints)))
      }
      l.gamma.phi<-function(parm){
        gamma=parm[1:(p[3])]; phi=parm[p[3]+1]
        return(mean(l.gamma.phi.each(gamma,phi)))
      }
      env.l.gamma.phi=new.env()
      assign("gamma", gamma, envir = env.l.gamma.phi)
      assign("phi", phi, envir = env.l.gamma.phi)
      d1l.gamma.phi<-attr(numericDeriv(quote(l.gamma.phi.each(gamma,phi)), c('gamma','phi'), env.l.gamma.phi),'gradient')

      d2l.gamma.phi<-hessian(l.gamma.phi,x=c(gamma,phi))
      colnames(d1l.gamma.phi)<-c(paste0(rep('gamma',p[3]),1:p[3]),'phi'); colnames(d2l.gamma.phi)=colnames(d1l.gamma.phi); rownames(d2l.gamma.phi)=colnames(d1l.gamma.phi)
    }
    if(verbose){print('1st Derivative of l(gamma,phi)')}
    if(verbose){print('2nd Derivative of l(gamma,phi)')}

    # calculate covariance matrix for beta
    gj<-aggregate(gij~id,FUN=sum,data=data.frame(id,gij))[,2:(p[1]+1)]
    inside<-solve(dsm.beta)%*%(t(gj)-d1.gij[,grepl('alpha',colnames(d1.gij))|colnames(d1.gij)=='tau']%*%solve(d2l.alpha.tau)%*%t(d1l.alpha.tau)-d1.gij[,grepl('gamma',colnames(d1.gij))|colnames(d1.gij)=='sigma'|colnames(d1.gij)=='phi']%*%solve(d2l.gamma.phi)%*%t(d1l.gamma.phi))
    var.beta<-inside%*%t(inside)
    colnames(var.beta)<-c(paste0(rep('beta',p[1]),0:(p[1]-1))); rownames(var.beta)=colnames(var.beta)
    return(var.beta)
  }
  if(se==TRUE){var.beta<-SE(beta,alpha,gamma,sigma,tau,phi)}
  if(se==FALSE){var.beta=NULL}
  fitvalues<-c(as.matrix(x)%*%as.matrix(beta,ncol=1))
  gammalist<-gammalist
  result=list(Method="lmeaipw",Call=arg_checks,nr.conv=(max(abs(diff.nr))<=conv),
              nr.iter=nr.iter,nr.diff=max(abs(diff.nr)),beta=beta,var.beta=var.beta,fitvalues=fitvalues,gammalist=gammalist)
  return(result)
}
