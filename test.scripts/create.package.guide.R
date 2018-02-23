# load libraries
library(devtools)
library(roxygen2)

# set working directory to root of the project:
# /interplay.gut.microbiome.immune.system/
#     /development/
#     /elastic.net.regularization/
#     /junk/
#     /least.absolute.shrinkage/selection.operator/
#     /ridge.regression/

# set the root directory
setwd("/interplay.gut.microbiome.immune.system")

# create folder for new package
create("test")

# write function belonging to particular function with documentation,
# please use template:

#' Add together two numbers.
#' 
#' @author Matthijs Knigge
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @export
#' @examples
#' add(1, 1)
#' add(10, 1)
add <- function(x, y) {
  x + y
}

# save function to package of interest under R, for example:
#       /elastic.net.regularization/R/add.R

# then process documentation, first set package folder:
setwd("elastic.net.regularization")

# build documentation with: 
devtools::document()


