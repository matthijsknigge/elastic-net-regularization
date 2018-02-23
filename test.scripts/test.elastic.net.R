# libraries
library(MASS)  # Package needed to generate correlated precictors
library(glmnet)  # Package to fit ridge/lasso/elastic net models
library(ggplot2)

# seed
set.seed(19875)  # Set seed for reproducibility

# globals
n <- 1000  # Number of observations
p <- 5000  # Number of predictors included in model
real_p <- 15  # Number of true predictors
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- apply(x[,1:real_p], 1, sum) + rnorm(n)

# split data 
train_rows <- sample(1:n, .66*n) # random data set
x.train <- x[train_rows, ] # train (2/3)
x.test <- x[-train_rows, ] # test (1/3)

y.train <- y[train_rows]
y.test <- y[-train_rows]

# fit models





