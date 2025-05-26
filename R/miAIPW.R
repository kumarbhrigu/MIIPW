#' miAIPW
#' @title
#' Fit a geeglm model using miAIPW
#' @description provides augmented inverse probability weighted estimates of parameters for GEE model
#' of response variable using different covariance structure. The augmented terms are estimated by using multiple imputation model.
#' @details It uses the augmented inverse probability weighted method to reduce the bias
#' due to missing values in GEE model for longitudinal data. The response variable \eqn{\mathbf{Y}} is related to the coariates as \eqn{g(\mu)=\mathbf{X}\beta}, where \code{g} is the link function for the glm. The estimating equation is
#' \deqn{\sum_{i=1}^{k}\sum_{j=1}^{n}(\frac{\delta_{ij}}{\pi_{ij}}S(Y_{ij},\mathbf{X}_{ij},\mathbf{X}'_{ij})+(1-\frac{\delta_{ij}}{\pi_{ij}})\phi(\mathbf{V}=\mathbf{v}))=0}
#' where \eqn{\delta_{ij}=1} if there is missing value in covariates and 0 otherwise,
#' \eqn{\mathbf{X}} is fully observed all subjects and \eqn{\mathbf{X}'} is partially missing,
#'  where \eqn{\mathbf{V}=(Y,\mathbf{X})}. The missing score function values due to incomplete data are estimated
#'  using an imputation model through mice which we have considered as \eqn{\phi(\mathbf{V}=\mathbf{v}))}. The estimated value \eqn{\phi(\mathbf{V}=\mathbf{v}))} is obtained
#'  through multiple imputation.
#'
#' @param data longitudinal data set where each subject's outcome has been measured at same time points and number
#'  of visits for each patient is similar.
#' Covariance structure of the outcome variable like "unstuctured","independent","AR1"
#' ,"Exchageable"
#' @param formula formula for the response model
#' @param id column name of id of subjects in the dataset
#' @param visit column name of timepoints of visit in the dataset
#' @param family name of the distribution for the response variable, For more information on how to use \code{family} objects, see \code{\link[stats]{family}}
#' @param init.beta initial values for the regression coefficient of GEE model
#' @param init.alpha initial values for the correlation structure
#' @param init.phi initial values for the csale parameter for
#' @param tol tolerance in calculation of coefficients
#' @param weights A vector of weights for each observation.
#' If an observation has weight 0, it is excluded from the calculations of any parameters. Observations with a NA anywhere (even in variables not included in the model) will be assigned a weight of 0. Weights are updated as the mentioned the details.
#' @param corstr a character string specifying the correlation structure. It could "independent", "exchangeable", "AR-1", "unstructured"
#' @param maxit maximum number iteration for newton-raphson
#' @param m number of imputation used to update the missing score function value due incomplete data.
#' @param pMat predictor matrix as obtained in \code{\link[mice]{mice}}
#' @param method method option for mice model,for information see \link[mice]{mice}
#' @return A list of objects containing the following objects
#' \describe{
#'   \item{call}{details about arguments passed in the function}
#'   \item{beta}{estimated regression coeffictient value for the response model}
#'   \item{niter}{number of iteration required}
#'   \item{betalist}{list of beta values at different iteration}
#'   \item{weight}{estimated weights for the observations}
#'   \item{mu}{mu values according \link[stats]{glm}}
#'   \item{phi}{etsimated phi value for the \code{glm} model}
#'   \item{hessian}{estimated hessian matrix obtained from the last iteration}
#'   \item{betaSand}{sandwich estimator value for the variance covariance matrix of the beta}
#' }
#' @import mice
#' @import MASS
#' @import Matrix
#' @export
#' @references Wang, C. Y., Shen-Ming Lee, and Edward C. Chao. "Numerical equivalence of imputing scores and weighted estimators in regression analysis with missing covariates." Biostatistics 8.2 (2007): 468-473.
#' @references Seaman, Shaun R., and Stijn Vansteelandt. "Introduction to double robust methods for incomplete data." Statistical science: a review journal of the Institute of Mathematical Statistics 33.2 (2018): 184.
#' @references Vansteelandt, Stijn, James Carpenter, and Michael G. Kenward. "Analysis of incomplete data using inverse probability weighting and doubly robust estimators." Methodology: European Journal of Research Methods for the Behavioral and Social Sciences 6.1 (2010): 37.
#' @examples
#'  \dontrun{
#' ##
#' formula<-C6kine~ActivinRIB+ActivinRIIA+ActivinRIIAB+Adiponectin+AgRP+ALCAM
#' pMat<-mice::make.predictorMatrix(srdata1[names(srdata1)%in%all.vars(formula)])
#' m1<-miAIPW(data=srdata1,
#' formula<-formula,id='ID',
#'  visit='Visit',family='gaussian',init.beta = NULL,
#' init.alpha=NULL,init.phi=1,tol=.00001,weights = NULL,
#' corstr = 'exchangeable',maxit=4,m=2,pMat=pMat)
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso   \link[MIIPW]{SIPW},\link[MIIPW]{miSIPW},\link[MIIPW]{miAIPW}
#'
miAIPW<-function(data,formula,id,visit,family,init.beta=NULL,init.alpha=NULL,
               init.phi=NULL,tol=0.001,weights=NULL,corstr='independent',maxit=50,m=2,pMat,method=NULL){
  call<-match.call()
  if(is.null(data)){
    stop("data can't be NULL")
  }
  if(is.null(id)){
    stop("id can't be NULL")
  }
  if(is.null(visit)){
    stop("visit can't be NULL")
  }
  if(is.null(family)){
    stop("family can't be NULL, must be from exponential distribution")
  }
  dt<-data
  stopifnot(id%in%names(dt))
  stopifnot(visit%in%names(dt))


  modeldata<-model.frame(formula,dt)
  names(dt)[which(names(dt)==id)]<-'id'
  names(dt)[which(names(dt)==visit)]<-'visit'

  dtid<-data.frame(table(dt$id))[2]
  dtidvisit<-table(dt$id,dt$visit)
  #if(max(dtid)>maxvisit){
   # stop("Some id's having visits more than",maxvisit)
  #}
  uniqueVisit<-unique(dt$visit)
  for(i in 1:length(uniqueVisit)){
    dt[dt$visit==uniqueVisit[i],]$visit<-i
  }
  init.beta<-init.beta;init.alpha<-init.alpha;init.phi<-init.phi

  # function to create IPW

  #mdata<-data.frame(id=c(rep(1,4),rep(2,4),rep(3,4)),visit=rep(1:4,3),y=rnorm(12,0,1),x1=rnorm(12,2,1),x3=c(1,1,1,1,1,NA,NA,NA,NA,NA,1,1))
  #data<-mdata;id<-'id';visit<-'visit';formula<-y~x1+x3
  #data<-dt;id<-'id';visit<-'visit';formula<-C6kine~ActivinRIB+ActivinRIIA+ActivinRIIAB+Adiponectin+AgRP+ALCAM
  #pMat<-make.predictorMatrix(data[names(data)%in%all.vars(formula)]);m=1;method=NULL
  if(anyNA(dt)==T){
   miaipwlong<-function(data,id,visit,formula)
  {
    dt<-data
    mterms<-all.vars(formula)
    #names(dt)[which(names(dt)==y)]<-'y'
    names(dt)[which(names(dt)==visit)]<-'visit'
    names(dt)[which(names(dt)==id)]<-'id'
    nvisit<-length(unique(dt$visit))
    k1<-which(colnames(dt)=='visit')
    k2<-which(colnames(dt)=='id')
    dt1<-dt[names(dt)%in%c('id','visit',mterms)==T]

    dt<-cbind(dt,r.ptrn=as.numeric(apply(dt,1,anyNA)))
    dt$r.ptrn<-ifelse(dt$r.ptrn==0,1,0)
    dtlist<-split(dt,dt$visit)
    mprob1<-NULL
    for(i in 1:nvisit){
      dt2<-dtlist[[i]]
      xptrn<-apply(dt2,2,anyNA)
      xmodel<-mterms
      fxmodel<-xptrn[names(xptrn)%in%mterms]
      fxmodel<-fxmodel[fxmodel==F]
      dt2new<-dt2[names(dt2)%in%mterms]
      micemodel<-mice::mice(
        data=dt2new,
        m=m,
        method=method,
        predictorMatrix=pMat,
        ignore=NULL,printFlag = F

      )
      completeList<-complete(micemodel,1)
      #names(dt2)[names(dt2)%in%names(fxmodel)]<-names(completeList)[names(dt2)%in%names(fxmodel)]
      dt2[names(dt2)%in%mterms]<-completeList
      modeldt2<-formula(paste0('r.ptrn~',paste(names(fxmodel),collapse = '+')))
      m1<-glm(modeldt2,data=dt2,family='binomial')
      dt2<-data.frame(dt2,w=m1$fitted.values)
      for(j in 1:length(dt2$id)){
        if(dt2$r.ptrn[j]==1){
          dt2$w[j]<-1/dt2$w[j]
        }else{
          dt2$w[j]<-1/(1-dt2$w[j])
        }
      }
      #dt2$w<-sapply(dt2$r.ptrn,dt2$w,1-dt$w)
      mprob1[[i]]<-dt2
    }
    reddt<-Reduce('rbind',mprob1)
    reddt<-reddt[order(reddt$id),]
    reddt
  }
  wlong<-miaipwlong(data=dt,id=id,visit=visit,formula=formula)
  #weights<-rep(1,nrow(data))
  wlong[all.vars(formula)]<-dt[all.vars(formula)]
  ##
  pMat1=make.predictorMatrix(wlong[names(wlong)%in%all.vars(formula)])
  wlongMice<-mice(
    data=wlong,
    m=m,
    method=NULL,
    predictorMatrix=pMat1,
    ignore=NULL,printFlag = F
  )
  #cwlongMice<-complete(wlongMice)
  imputeId<-list()
  for(i in 1:m){
    imputeIdm<-complete(wlongMice,action=i)

    imputeId[[i]]<-imputeIdm[imputeIdm$r.ptrn==0,]
  }

  nIdm<-nrow(imputeId[[1]])
  imputeId<-Reduce('rbind',imputeId)
  idmat<-matrix(imputeId$id,nrow=nIdm)
  seqId<-seq(max(unique(wlong$id))+1,max(unique(wlong$id))+(ncol(idmat)-1)*length(unique(idmat[,1])))
  unFisrt<-rep(data.frame(table(idmat[,1]))$Freq,(ncol(idmat)-1))
  repseqId<-list()
  for(l in 1:length(seqId)){
    repseqId[[l]]<-rep(seqId[l],unFisrt[l])
  }
  repseqId<-unlist(repseqId)

  imputeId$id[-c(1:nIdm)]<-repseqId

  imputeId$w<-imputeId$w/m
  wlong<-rbind(wlong,imputeId)
  wlong$w<-ifelse(apply(wlong,1,anyNA),0,wlong$w)
  if(sum(wlong$r.ptrn)!=nrow(wlong)){
    weights<-wlong$w
  }else{
    weights<-rep(1,nrow(data))
  }

  dt<-na.omit(wlong)
  weights<-dt$w
  }else{
    dt<-dt
    weights<-rep(1,nrow(data))
  }

  dat <- model.frame(formula=formula, data=dt, na.action=na.exclude)
  len<-data.frame(table(dt$id))
  nn <- dim(dat)[1]
  corlist<-c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", "fixed", "userdefined")
  cor.match <- charmatch(corstr, corlist)

  fmly<-get(family, mode = "function", envir = parent.frame(2))
  fmly<-fmly()
  LinkFun <- fmly$linkfun
  InvLink <- fmly$linkinv
  VarFun <- fmly$variance
  InvLinkDeriv <- fmly$mu.eta

  modterms <- terms(formula)
  X <- model.matrix(formula,dat)
  Y <- model.response(dat)
  p <- dim(X)[2]
  W <- Diagonal(x=weights)
  sqrtW <- sqrt(W)
  K <- length(unique(dt$id))
  StdErr <- Diagonal(nn)
  dInvLinkdEta <- Diagonal(nn)
  Resid <- Diagonal(nn)


  if(is.null(init.phi)){
    init.phi<-1
  }else{
    init.phi<-init.phi
  }
  phi<-init.phi

  linkOfMean <- LinkFun(mean(Y))
  if(is.null(init.beta)){
    init.beta <- rep(0, dim(X)[2])
    init.beta[1] <- linkOfMean
  }else{
    init.beta<-init.beta
  }
  beta <- init.beta

  if(is.null(init.alpha)){
    init.alpha<-init.alpha
  }else{
    init.alpha<-init.alpha
  }

  stop <- F
  converged <- F
  count <- 0
  unstable <- F
  phiold <- phi
  kbeta <- beta
  alpha<-init.alpha
  R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
  betalist<-list()
  p=0
  while(!stop){
    p <- p+1
    off<-0
    eta <- as.vector(X %*% beta) + off

    mu <- InvLink(eta)
    diag(dInvLinkdEta)<-InvLinkDeriv(eta)
    diag(StdErr) <- sqrt(1/VarFun(mu))

    #R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
    #R.alpha.inv<-as.matrix(R.alpha.inv)
    UpdatePhi<-function(y,x,vfun,mu,w){
      ymu<-(as.matrix(y-mu,ncol=1,nrow=length(y),byrow=T))
      res<-vfun%*%(w%*%(ymu))
      Newphi<-sum(res*res)/(length(y)-ncol(x))
      Newphi
    }
    phi.new<-UpdatePhi(y=Y,x=X,vfun=StdErr,mu=mu,w=W)
    phi<-phi.new

    #y=Y;x=X;vfun=StdErr;mu=mu;w=W;phi=1;corstr = corstr;ni=len$Freq;mv=NULL;id=id$Pig;visit=visit$Time
    #y=Y;x=X;vfun=StdErr;mu=mu;w=W;phi=phi;corstr = corstr;ni=len$Freq;mv=NULL;id=dt$id;visit=dt$visit

    updateALpha<-function(y,x,vfun,mu,w,phi,corstr,ni,mv=NULL,id,visit){
      nk<-max(ni)
      Alpha<-NULL
      Resid<-vfun%*%((y-mu))
      dt2<-cbind.data.frame(Res=as.vector(Resid),id=id,visit=visit)
      Ralpha<-matrix(1,ncol=nk,nrow=nk)
      ResidMatlist<-list()
      length(ResidMatlist)<-length(ni)
      for(i in 1:length(ResidMatlist)){
        a<-rep(NA,nk)
        dt1<-dt2[dt2$id==unique(dt2$id)[1],]
        a[c(which(c(1:nk)%in%dt1$visit==T))]<-dt1$Res
        a[is.na(a)]<-0
        ResidMatlist[[i]]<-triu(a%*%t(a))
      }
      resMAtsum<-Reduce('+',ResidMatlist)

      if(corstr=='fixed'){
        Ralpha<-Ralpha.fixed
      }else if(corstr=="independent"){
        for(i in 1:nk){
          for(j in 1:nk){
            if(i!=j){
              Ralpha[i,j]<-0
            }else{
              Ralpha[i,j]<-1
            }
          }
        }

      }else if(corstr=='unstructured'){
        Ralpha<-resMAtsum/((sum(ni)-ncol(x))*phi)
        diag(Ralpha)<-1
        #ResidNa<-bdiag(ResidMatlist)
        #ResidNa2<-band(ResidNa,k1=1,k2=dim(ResidNa)[2])

      }else if(corstr=='exchangeable'){
        Alpha<-(sum(resMAtsum)-sum(diag(resMAtsum)))/(phi*(sum(ni*ni)-sum(ni)-ncol(x)))
        for(i in 1:nk){
          for(j in 1:nk){
            if(i!=j){
              Ralpha[i,j]<-Alpha
            }else{
              Ralpha[i,j]<-1
            }
          }
        }
      }else if(corstr=='AR-1'){
        Alpha<-(sum(triu(resMAtsum,k=1)))/((sum(ni-1)-ncol(x))*phi)
        for(i in 1:nk){
          for(j in 1:nk){
            if(i!=j){
              Ralpha[i,j]<-Alpha^(abs(i-j))
            }else{
              Ralpha[i,j]<-1
            }
          }
        }
      }
      Ralpha
    }
    ralpha<-updateALpha(y=Y,x=X,vfun=StdErr,mu=mu,w=W,phi=phi,corstr = corstr,ni=len$Freq,mv=NULL,id=dt$id,visit=dt$visit)
    ralpha<-ralpha*phi
    ralpha.list<-list()
    blockdiag<-c()
    for(i in 1:length(unique(dt$id))){

      dt1<-dt[dt$id==unique(dt$id)[i],]
      ralpha.list[[i]]<-ginv(as.matrix(ralpha[c(dt1$visit),c(dt1$visit)]))
      blockdiag[i]<-length(dt1$id)
    }
    R.alpha.inv<-(bdiag(ralpha.list))
    #y=Y;x=X;vfun=StdErr;mu=mu;w=W;D=dInvLinkdEta;Ralpha=R.alpha.inv;beta=beta
    updateBeta<-function(y,x,vfun,mu,w,D,Ralpha,beta){
      hess<-t(D%*%x)%*%vfun%*%(Ralpha)%*%vfun%*%D%*%w%*%x
      u<-t(D%*%x)%*%vfun%*%(Ralpha)%*%vfun%*%w%*%as.matrix(y-mu)
      Newbeta<-beta+solve(as.matrix(hess),u)
      rslt<-list()
      rslt$Newbeta<-Newbeta
      rslt$hess<-hess
      rslt
    }
    Beta<-updateBeta(y=Y,x=X,vfun=StdErr,mu=mu,w=W,D=dInvLinkdEta,Ralpha=R.alpha.inv,beta=beta)
    betaNew<-Beta$Newbeta
    phi <- phi.new
    R.alpha.inv<-as.matrix(R.alpha.inv)
    phiold <- phi
    beta<-betaNew
    if( max(abs((beta - kbeta)/(beta+.Machine$double.eps ))) < tol ){converged <- T; stop <- T}
    if(p >= maxit){stop <- T}

    betalist[[p]]<-betaNew
    betaold<-beta
  }
  eta <- as.vector(X %*% beta)
  mu <- InvLink(eta)
  diag(dInvLinkdEta)<-InvLinkDeriv(eta)
  diag(StdErr) <- sqrt(1/VarFun(mu))
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

    sigmahat<-Blockdiag
    usand<-t(x)%*%t(D)%*%vfun%*%Ralpha%*%vfun%*%w%*%sigmahat%*%w%*%vfun%*%Ralpha%*%vfun%*%D%*%x
    hessmat<-as.matrix(hessmat)
    upSandwich<-t(solve(hessmat,usand))
    #upSandwich<-upSandwich%*%solve(hessmat)
    upSandwich<-t(solve(t(hessmat),upSandwich))
    upSandwich
  }
  betaSand<-updateSandW(y=Y,x=X,vfun=StdErr,mu=mu,w=W,D=dInvLinkdEta,Ralpha=R.alpha.inv,beta=beta,
                        hessmat = Beta$hess,blockdiag=blockdiag)
  reslt<-list()
  reslt$call<-call
  reslt$beta<-beta
  reslt$niter<-p
  reslt$betalist<-betalist
  reslt$Ralpha<-ralpha/phi
  reslt$weight<-W
  reslt$mu<-mu
  reslt$y<-Y
  reslt$phi<-phi
  reslt$hessian<-Beta$hess
  reslt$betaSand<-betaSand
  class(reslt) <- "ipw"
  return(reslt)


}

utils::globalVariables(c("Ralpha.fixed","ginv","triu","model.frame", "model.matrix",
                         "model.response", "na.exclude", "na.omit", "na.pass",
                         "pnorm", "terms","stats","glm"))




