#' Inverse Rank
#' @author Matthijs Knigge
#'
#' @param phenotype numeric vector
#' @keywords plot
#' @export
#' @examples
#' glmnet.inverse.rank()
#' @return forces normal distribution of vector
#' 
#'      
glmnet.inverse.rank <- function(phenotype){
  # libraries
  
  # residuals
  res <- rank(phenotype)
  # transform to normal distribution
  res <- qnorm(res/(length(res)+0.5))
  # return vector
  return(res)
}
