
#' summary method for ipw
#'
#' @param object ipw object
#' @param ... further argument can be passed
#'
#' @return summary of ipw object
#'
#' @export
summary_ipw<- function(object, ...)  {
  Coef<-list()
  Coef[[1]] <- c(as.vector(object$beta))
  smat<-as.matrix(object$betaSand)
  colnames(smat)<-NULL;row.names(smat)<-NULL
  Coef[[2]] <- sqrt(diag(smat))
  Coef[[3]] <- as.vector(Coef[[1]])/as.vector(Coef[[2]])
  Coef[[4]] <- round(2*pnorm(abs(Coef[[3]]), lower.tail=F), digits=8)
  dtCoef<-data.frame(Coef[[1]],Coef[[2]],Coef[[3]],Coef[[4]])
  colnames(dtCoef) <- c("Estimates","SE", "z value", "Pr(>|z|)")
  coefnames<-rownames(object$beta)
  summ <- list( call=object$call,inference=dtCoef,
                phi = object$phi,Ralpha=object$Ralpha)
  class(summ) <- 'summary_ipw'
  return(summ)
}

#' summary method for meanscore
#'
#' @param object meanscore object
#' @param ... further argument can be passed
#'
#' @return summary of meanscore object
#' @export
summary_meanscore<- function(object, ...)  {
  Coef<-list()
  Coef[[1]] <- c(as.vector(object$beta))
  smat<-as.matrix(object$betaSand)
  colnames(smat)<-NULL;row.names(smat)<-NULL
  Coef[[2]] <- sqrt(diag(smat))
  Coef[[3]] <- as.vector(Coef[[1]])/as.vector(Coef[[2]])
  Coef[[4]] <- round(2*pnorm(abs(Coef[[3]]), lower.tail=F), digits=8)
  dtCoef<-data.frame(Coef[[1]],Coef[[2]],Coef[[3]],Coef[[4]])
  colnames(dtCoef) <- c("Estimates","SE", "z value", "Pr(>|z|)")
  coefnames<-rownames(object$beta)
  summ <- list( call=object$call,inference=dtCoef,
                phi = object$phi,Ralpha=object$Ralpha)
  class(summ) <- 'summary_meanscore'
  return(summ)
}

