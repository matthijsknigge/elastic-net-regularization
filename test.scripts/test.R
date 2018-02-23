# libraries

library(ggplot2)
library(cross.validation)
library(elastic.net.regularization)


# Generate data
set.seed(19874)
n <- 10    # Number of observations
p <- 500     # Number of predictors included in model
real_p <- 150  # Number of true predictors
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- apply(x[,1:real_p], 1, sum) + rnorm(n)
nrow(x)
# add colnames
rownames(x) <- sprintf("sample%s", seq(1:nrow(x)))
# vector names
y <- setNames(y, sprintf("sample%s", seq(1:length(y))))
# --------------------------------------------


result <- cv.glmnet.wrapper(y = y, x = x, kfold = 3)
