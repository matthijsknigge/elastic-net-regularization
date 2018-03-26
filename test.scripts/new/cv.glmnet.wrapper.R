#' Add together two numbers.
#'
#' @author Olivier
#' @param x x matrix
#' @param y a vector of outcomes.
#' @param k an integer for the number of folds.
#' @return Elastic Net Results
#' @export
#' @examples
#' cv.glmnet.wrapper()
cv.glmnet.wrapper <- function(y, x, kfold=10, scale.x=F, verbose=T,  ...) {
  # Scale & center the predictor set
  if (scale.x) {
    x           <- x[,apply(x, 2, var) != 0]
    x           <- scale(x)
  }

  # Create the 10 stratified folds
  folds         <- create.k.folds(y, k=kfold)
  # Run glmnet.wrapper on the k folds
  results <- lapply(folds, function(fold, y, x) {
    # apply glmnet for every fold
    return(glmnet.wrapper(y=y[-fold],
                          x=x[-fold, , drop=F],
                          newx=x[fold,  , drop=F],
                          newy=y[fold],
                          ...))
  }, y=y, x=x)
  return(results)
}



