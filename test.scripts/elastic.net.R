## ------------------------------------------------------------------------
## Wrapper for glmnet which has some added features and automaticly
## chooses the best alpha & lambda combination.
glmnet.wrapper <- function(y, x, newx=NULL, newy=NULL, alpha.runs=NULL,
                           parallel=F, return.type="custom",
                           method="spearman", use="pairwise.complete.obs",
                           remove.outliers=F) {
  # Remove NA values
  y             <- y[!is.na(y)]
  if (remove.outliers) {
    # Remove outliers in the response
    # (improves prediction but reduces sample size)
    y <- y[!y %in% boxplot.stats(y)$out]
  }
  ol            <- intersect(rownames(x), names(y))
  y             <- y[ol]
  x             <- as.matrix(x[ol, , drop=F])
  models.alpha  <- list()
  alpha.vec     <- c()
  if (is.null(alpha.runs)) {
    alpha.runs    <- c(0, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.65, 0.9, 1)
  }
  # Loop over the possible alpha values and fit the model
  for (a in 1:length(alpha.runs)) {
    model             <- cv.glmnet(y=y, x=x, alpha=alpha.runs[a], parallel=parallel)
    #alpha.vec[a]      <- model$cvm[model$lambda == model$lambda.min]
    # Use 1se to prevent overfitting
    alpha.vec[a]      <- model$cvm[model$lambda == model$lambda.1se]
    models.alpha[[a]] <- model
  }
  # Arrange the results in a custom list format which contains predictions
  if (return.type=="custom") {
    result                <- list()
    best.alpha            <- which(alpha.vec == min(alpha.vec))
    best.model            <- models.alpha[[best.alpha]]
    best.lambda           <- best.model$lambda.1se
    coef                  <- as.matrix(coef(best.model, s=best.lambda))
    coef                  <- coef[coef != 0,]
    if (!is.null(newx) && !is.null(newy)) {
      pred                <- predict(best.model$glmnet.fit, newx=as.matrix(newx),
                                     s=best.lambda, type="response")
      pred.names          <- rownames(pred)
      pred                <- as.numeric(pred)
      names(pred)         <- pred.names
      pred.cor            <- cor(pred, newy, method=method, use=use)
      result$pred         <- pred
      result$pred.cor     <- pred.cor
      result$n.test       <- length(newy)
      result$mse          <- sum((pred - newy)^2)/length(pred)
    }
    result$best.alpha      <- alpha.runs[best.alpha]
    result$best.lambda     <- best.lambda
    result$best.model      <- best.model
    result$coef            <- coef
    result$n.train         <- length(y)
  } else if (return.type == "model") {
    result <- models.alpha[[which(alpha.vec == min(alpha.vec))]]
  } else if (return.type == "all") {
    result <- models.alpha
  }
  return(result)
}
## ------------------------------------------------------------------------
## K-fold cross validation for glmnet.wrapper.
cv.glmnet.wrapper <- function(y, x, kfold=10, scale.x=F, verbose=T,  ...) {
  # Scale & center the predictor set
  if (scale.x) {
    x           <- x[,apply(x, 2, var) != 0]
    x           <- scale(x)
  }
  # Remove NA values
  y             <- y[!is.na(y)]
  x             <- na.omit(x)
  ol            <- intersect(rownames(x), names(y))
  y             <- y[ol]
  x             <- x[ol, , drop=F]
  # Create the 10 stratified folds
  folds         <- createFolds(y, k=kfold)
  # Run glmnet.wrapper on the k folds
  results <- lapply(folds, function(fold, y, x) {
    # This was very memory inefficient
    #x.train <- x[-fold, , drop=F]
    #x.test  <- x[fold,  , drop=F]
    #y.train <- y[-fold]
    #y.test  <- y[fold]
    return(glmnet.wrapper(y=y[-fold],
                          x=x[-fold, , drop=F],
                          newx=x[fold,  , drop=F],
                          newy=y[fold],
                          ...))
  }, y=y, x=x)
  return(results)
}