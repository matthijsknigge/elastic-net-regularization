#' Predict immune ohenotypes
#' @author Matthijs Knigge
#'
#' @param x matrix n*m, must have row, and colnames.
#' @param coef matrix n*m, must have row, and colnames.
#' @keywords predict
#' @export
#' @examples
#' predict.coef()
#' @return prediceted immune phenotypes
#' 
#'  
predict.coef <- function(x, coef) {
  # test for intersect
  if (length(intersect(names(coef), colnames(x))) != (length(coef))) {
    # intersect
    int <- intersect(names(coef), colnames(x))
    # apply intersect x
    x <- as.matrix(cbind(1, x[,names(coef[int])]))
    # multiply
    return(as.vector(coef[c("(Intercept)",int)]) %*% t(x))
  }
  # if statement succes
  else{
    # intersect
    int <- names(coef)[-1]
    # cbind data
    x <- as.matrix(cbind(1, x[,names(coef[int])]))
    # return
    return(as.vector(coef[c("(Intercept)",int)]) %*% t(x))
  }
}