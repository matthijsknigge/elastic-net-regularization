#' Parse elastic net results
#' @author Matthijs Knigge
#'
#' @param dir character, absolute path to root folder of results
#' @param data.y matrix n*m, must have row, and colnames. Phenotype data
#' @param folds numeric, amount of folds used
#' @param merge boolean, add group column to frame
#' @param by matrix 2*n, character. containing identifer and group column
#' @param melt boolean, default FALSE. melt output
#' @param count.input boolean, default FALSE. count amount of genes and SNPs that are input.
#' @param count.ouput boolean, default FALSE. count amount of genes and SNPs that are used.
#' @keywords plot
#' @export
#' @examples
#' glmnet.parse.results()
#' @return results
#'
#'
glmnet.parse.results <- function(dir, folds, data.y, merge = FALSE, data.c, melt = FALSE, count.input = FALSE, count.output = FALSE){
  # libraries
  require(data.table)

  # define data.frame for correlations
  df.coef <- data.frame(imputation = character(), identifier = character(), stimulation = character(),
                        coef.fold.1 = numeric(), coef.fold.2 = numeric(), coef.fold.3 = numeric(),
                        coef.fold.4 = numeric(), coef.fold.5 = numeric(), coef.fold.6 = numeric(),
                        coef.fold.7 = numeric(), coef.fold.8 = numeric(), coef.fold.9 = numeric(),
                        coef.fold.10 = numeric(), genes.input = numeric(), genes.output = numeric(),
                        snps.input = numeric(), snps.output = numeric(), alpha.1 = numeric(), alpha.2 = numeric(), 
                        alpha.3 = numeric(), alpha.4 = numeric(), alpha.5 = numeric(), alpha.6 = numeric(),
                        alpha.7 = numeric(), alpha.8 = numeric(), alpha.9 = numeric(), alpha.10 = numeric(), coef.mean = numeric())

  # define data.frame for pred ~ obs
  df.pred.obs <- data.frame(pred = numeric(), obs = numeric(), fold = character(),
                            imputation = character(), identifier = character(),
                            gene.amount = numeric())

  # list all phenotypes
  phenotypes <- list.dirs(dir, full.names = F, recursive = F)

  # for every phenotype
  for(pheno in phenotypes){
    # determine y
    y <- data.y[[pheno]]; names(y) <- rownames(data.y)
    # for fold in cytokine
    for(i in 1:folds){
      # if model exists
      if(dir.exists(file.path(paste0(dir, pheno, "/", "fold.", i, "/")))){
        # if file exists
        if (file.exists(paste0(dir, pheno, "/", "fold.", i, "/", pheno, ".rda"))) {
          # load elastic net
          load(paste0(dir, pheno, "/", "fold.", i, "/", pheno, ".rda"))
          # vector for coeff
          coeff.vector <<- c()
          genes.input  <<- c()
          genes.output <<- c()
          snps.input  <<- c()
          snps.output <<- c()
          alpha.vec <<- c()
          # for models in elastic net
          for(m in 1:length(result)){
            # get model from elastic net
            model <- result[[m]]
            # get spearman coefficient
            spearman.corr <- model$pred.cor
            # get predicted
            pred <- model$pred
            # get observed, need gene expression data
            obs  <- y[names(model$pred)]
            # get cross number
            imputation <- paste0("imputation.", i)
            # get fold nuber
            fold  <- paste0("fold.", m)
            # get cytokine name
            pheno.name <- pheno
            
            # amount SNPs as input
            snps.input  <<- c(snps.input, model$genes.amount)
            # amount of genes as input
            genes.input <<- c(genes.input, model$snps.amount)
            # amount of SNPs used in model
            snps.output <<- c(snps.output, length(names(model$coef)[grepl(pattern = "rs", x = names(model$coef))]))
            # amount of genes used in model
            genes.output<<- c(genes.output, length(names(model$coef)[!grepl(pattern = "rs", x = names(model$coef))]))
            # bind spearman correlation
            coeff.vector <<- c(coeff.vector, spearman.corr)
            # bind alpha
            alpha.vec <<- c(alpha.vec, model$best.alpha)
            # bind info together obs and pred
            df.pred.obs <- rbind(df.pred.obs, data.frame("pred" = pred,
                                                         "obs" = obs,
                                                         "fold" =fold,
                                                         "imputation" = imputation,
                                                         "identifier" = pheno.name))

          }
          df.coef <- rbind(df.coef, data.frame("imputation" = imputation, "identifier" = pheno.name,
                                               "coef.fold.1" = coeff.vector[1], "coef.fold.2" = coeff.vector[2], "coef.fold.3" = coeff.vector[3], "coef.fold.4" = coeff.vector[4],
                                               "coef.fold.5" = coeff.vector[5], "coef.fold.6" = coeff.vector[6], "coef.fold.7" = coeff.vector[7], "coef.fold.8" = coeff.vector[8],
                                               "coef.fold.9" = coeff.vector[9], "coef.fold.10" = coeff.vector[10], 
                                               "genes.input" = mean(genes.input, na.rm = T), 
                                               "genes.output" = mean(genes.output, na.rm = T),
                                               "snps.input" = mean(snps.input, na.rm = T),
                                               "snps.output" = mean(snps.output, na.rm = T),
                                               "alpha.1" = alpha.vec[1], "alpha.2" = alpha.vec[2], "alpha.3" = alpha.vec[3], "alpha.4" = alpha.vec[4],
                                               "alpha.5" = alpha.vec[5], "alpha.6" = alpha.vec[6], "alpha.7" = alpha.vec[7], "alpha.8" = alpha.vec[8],
                                               "alpha.9" = alpha.vec[9], "alpha.10" = alpha.vec[10], coef.mean = mean(coeff.vector)))
          alpha.vec <<- NULL
          coeff.vector <<- NULL
          genes.input  <<- NULL
          genes.output <<- NULL
          snps.input  <<- NULL
          snps.output <<- NULL
        }
      }
    }
  }

  # remove from global environment
  rm(coeff.vector, envir = .GlobalEnv); rm(genes.input, envir = .GlobalEnv); rm(genes.output, envir = .GlobalEnv); rm(snps.input, envir = .GlobalEnv); rm(snps.output, envir = .GlobalEnv);

  # merge results
  if(merge){
    df.pred.obs <- merge(df.pred.obs, data.c, by = "identifier")
    df.coef <- merge(df.coef, data.c, by = "identifier")

    # melt coef
    if(melt){
      # create new frame
      df.coef.new <- data.frame(identifier = character(), group = character(), coef = numeric())

      # for id in coefficients
      for(id in unique(df.coef$identifier)){
        # subset
        tmp <- df.coef[which(df.coef$identifier == id), ]
        # create vector from coefficients
        tmp.vec   <- as.vector(as.matrix(tmp[,c("coef.fold.1", "coef.fold.2", "coef.fold.3", "coef.fold.4", "coef.fold.5", "coef.fold.6","coef.fold.6", "coef.fold.7", "coef.fold.8", "coef.fold.9", "coef.fold.10")]))
        # create alpha vector
        tmp.alpha <- as.vector(as.matrix(tmp[, c("alpha.1", "alpha.2", "alpha.3", "alpha.4", "alpha.5", "alpha.6","alpha.6", "alpha.7", "alpha.8", "alpha.9", "alpha.10")]))
        # bind to frame
        df.coef.new <- rbind(df.coef.new, data.frame("identifier" = id, "group" = unique(tmp$group),
                                               "coef" = tmp.vec, "alpha" = tmp.alpha))
      }
    }
    # bind to frame
    df.coef.melt <- df.coef.new
  }



  # return results
  return(list(df.coef = df.coef, df.pred.obs = df.pred.obs, df.coef.melt = df.coef.melt))

}
