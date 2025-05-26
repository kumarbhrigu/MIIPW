
#' print method for ipw
#'
#' @param x ipw object
#' @param ... further argument can be passed
#'
#' @return print result for ipw object
#' @export
print_ipw <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(x$beta))
  cat("\n Scale Parameter: ", signif(x$phi, digits=4), "\n")
  #if(x$corstr=="unstructured"){
  cat("\n Estimated Correlation: ","\n")
  print(x$Ralpha,digits=4)
  #}
  #else cat("\n Estimated Correlation: ", signif(x$rho, digits=4), "\n")

}


#' print method for meanscore
#'
#' @param x meanscore object
#' @param ... further argument can be passed
#'
#' @return print result for meanscore object
#' @export
print_meanscore <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(x$beta))
  cat("\n Scale Parameter: ", signif(x$phi, digits=4), "\n")
  #if(x$corstr=="unstructured"){
  cat("\n Estimated Correlation: ","\n")
  print(x$Ralpha,digits=4)
  #}
  #else cat("\n Estimated Correlation: ", signif(x$rho, digits=4), "\n")

}
