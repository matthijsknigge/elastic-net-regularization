#' Simple Linear Regression
#' @author Matthijs Knigge
#'
#' This is a test function for personal purpose only, reflecting 
#' template that should be used. It sole purpose is just a reminder.
#' Simple linear regression: Y = B0 + B1 * X1 + E
#' @param x vector of length n, with independent variables. Default is cars data set speed.
#' @param y vector of length n, with dependent variables. Default is cars data set dist.
#' @param data n*p matrix with predixtor variables. Default is cars data set.
#' @keywords linear regression
#' @export
#' @examples
#' simple.linear.regression(x = c(1,2,3,4,5,6), y = c(1,2,3,4,5,6))
#' @return List with the following elements:
#'         height: coefficient
#'         intercept: coefficient
#'         height.se: unit increase in X, expect that Y increases by height.se
#'         intercept.se: unit increase in Y, expect that X increases by intercept.se
#'         height.p: significance of change
#'         intercept.p: significance of change
#'         
simple.linear.regression <- function(y = cars$dist, x = cars$speed, data = cars){
  # perform linear regression
  fit <- lm(y ~ x, data = data)
  # catch summary
  smmry <- summary(fit)
  # height
  height <- smmry$coefficients[2]
  # intercept
  intercept <- smmry$coefficients[1]
  # height se
  height.se <- smmry$coefficients[2, 2]
  # intercept se
  intercept.se <- smmry$coefficients[1, 2]
  # height p-value
  height.p <- smmry$coefficients[2, 4]
  # intercept p-value
  intercept.p <- smmry$coefficients[1, 4]
  # return parameters
  return(list(height = height,
              intercept = intercept,
              height.se = height.se,
              intercept.se = intercept.se,
              height.p = height.p,
              intercept.p = intercept.p))
}
