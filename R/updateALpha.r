#' @title internal function for updating alpha
#' @details arguments are from Fisher Scoring Algorithm
#' @param y response value for GEE model
#' @param x model matrix for the GEE model
#' @param vfun variance function for the GLM
#' @param mu mu vector for the GLM
#' @param w weight matrix
#' @param phi scale parameter
#' @param corstr correlation structure
#' @param ni list of visits per subject
#' @param mv NULL
#' @param id id column
#' @param visit visit column
#'
#' @export
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
          Ralpha[i,j]<-Alpha^(j-1)
        }else{
          Ralpha[i,j]<-1
        }
      }
    }
  }
  Ralpha
}

utils::globalVariables(c("Ralpha.fixed","ginv","triu"))
