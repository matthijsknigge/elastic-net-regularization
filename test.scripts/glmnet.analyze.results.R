# libraries
library(data.table)

# globals
folds <- 100
cytokine.levels <- "/groups/umcg-wijmenga/tmp04/data/500FG/cytokine_levels/500FG_log2_cytokine_levels.txt"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/500FG_RNAseq_TMM_normalized_read_counts.txt"
path.to.results <-"/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/"

# globals for testing
cytokine.levels <- "~/Bioinformatics/data/cytokine.levels.test.txt"
folds <- 10
cytokine.levels.info <- "~/Bioinformatics/data/cytokine.levels.info.txt"
gene.expression <- "~/Bioinformatics/data/expression.test.txt"
path.to.results <- "~/Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/results/"

# COLLECT DATA

# cytokine levels
c.l   <- fread(cytokine.levels); rownames(c.l) <- c.l$id; c.l$id <- NULL
# gene expression
g.e   <- fread(gene.expression, data.table = T); rownames(g.e) <- g.e$id; g.e$id <- NULL; g.e <- data.matrix(g.e)
# x
x <- g.e
# read cytokine info
y.info <- fread(cytokine.levels.info)

# define data.frame for correlations
df.coef <- data.frame(imputation = character(), cytokine = character(), stimulation = character(), coef.fold.1 = numeric(), coef.fold.2 = numeric(), coef.fold.3 = numeric(), coef.fold.4 = numeric(), 
                      coef.fold.5 = numeric(), coef.fold.6 = numeric(), coef.fold.7 = numeric(), coef.fold.8 = numeric(), coef.fold.9 = numeric(), 
                      coef.fold.10 = numeric())
# define data.frame for pred ~ obs
df.pred.obs <- data.frame(pred = numeric(), obs = numeric(), 
                          fold = character(), imputation = character(),
                          cytokine = character(), stimulation = character(), gene.amount = numeric())
# for every cytokine
for(cytokine in y.info$cytokine){
  # determine y
  y <- c.l[[cytokine]]; names(y) <- rownames(c.l)
  # for fold in cytokine
  for(i in 1:folds){
    # if model exists
    if(dir.exists(file.path(paste0(path.to.results, cytokine, "/", "fold.", i, "/")))){
      # if file exists
      if (file.exists(paste0(path.to.results, cytokine, "/", "fold.", i, "/", cytokine, ".rda"))) {
        # load elastic net
        load(paste0(path.to.results, cytokine, "/", "fold.", i, "/", cytokine, ".rda"))
        # vector for coeff
        coeff.vector <<- c()
        gene.amount  <<- c()
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
          cytokine.name <- cytokine
          # get stimulation
          stimulation <- y.info[which(y.info$cytokine == cytokine.name), ]$stimulation
          # amount of genes used
          x <- model$coef

          gene.amount <<- c(gene.amount, length(x[x > 0]))
          
          # bind spearman correlation
          coeff.vector <<- c(coeff.vector, spearman.corr)
          
          # bind info together obs and pred
          df.pred.obs <- rbind(df.pred.obs, data.frame("pred" = pred, 
                                                       "obs" = obs, 
                                                       "fold" =fold, 
                                                       "imputation" = imputation, 
                                                       "cytokine" = cytokine.name, 
                                                       "stimulation" = stimulation))
          
        }
        df.coef <- rbind(df.coef, data.frame("imputation" = imputation, "cytokine" = cytokine, "stimulation" = stimulation,
                                             "coef.fold.1" = coeff.vector[1], "coef.fold.2" = coeff.vector[2], "coef.fold.3" = coeff.vector[3], "coef.fold.4" = coeff.vector[4],
                                             "coef.fold.5" = coeff.vector[5], "coef.fold.6" = coeff.vector[6], "coef.fold.7" = coeff.vector[7], "coef.fold.8" = coeff.vector[8],
                                             "coef.fold.9" = coeff.vector[9], "coef.fold.10" = coeff.vector[10], "gene.amount" = mean(gene.amount)))
        coeff.vector <<- c()
        gene.amount  <<- c()
      }
    }
  }
}

df.pred.obs$samples <- rownames(df.pred.obs)
write.table(x = df.pred.obs, file = "Bioinformatics/test.results/data.summary/df.pred.obs.txt", quote = T, col.names = T, row.names = F)


# create new frame
coef.new <- data.frame(cytokine = character(), stimulation = character(), coef = numeric())

# for cytokine in coefficients
for(c in unique(df.coef$cytokine)){
  # subset
  tmp <- df.coef[which(df.coef$cytokine == c), ]
  # create vector from coefficients
  tmp.vec <- as.vector(as.matrix(tmp[,c("coef.fold.1", "coef.fold.2", "coef.fold.3", "coef.fold.4", "coef.fold.5", "coef.fold.6","coef.fold.6", "coef.fold.7", "coef.fold.8", "coef.fold.9", "coef.fold.10")]))
  # bind to frame
  coef.new <- rbind(coef.new, data.frame("cytokine" = c, "stimulation" = unique(tmp$stimulation),
                                         "coef" = tmp.vec))
}

write.table(x = coef.new, file = "Bioinformatics/test.results/data.summary/df.coef.molten.txt", quote = T, col.names = T, row.names = F)


new.df <- data.frame(cytokine = character(), stimulation = character(), amount.genes = numeric())
for(c in unique(df.coef$cytokine)){
  
  tmp <- df.coef[which(df.coef$cytokine == c), ]
  m <- mean(tmp$gene.amount)
  new.df <- rbind(new.df, data.frame(cytokine = unique(tmp$cytokine),
                                     stimulation = unique(tmp$stimulation),
                                     gene.amount = m))
  
}


write.table(x = new.df, file = "Bioinformatics/test.results/data.summary/cytokines.mean.genes.txt", quote = T, row.names = F, col.names = )




write.table(x = coef.new, file = "~/Bioinformatics/test.results/df.coef.new.txt", row.names = F, quote = T, col.names = T)


write.table(x = df.pred.obs, file = "~/Bioinformatics/test.results/df.pred.obs.txt", quote = T, row.names = T)
write.table(x = df.coef, file = "~/Bioinformatics/test.results/df.coef.txt", quote = T, row.names = F)


load("Bioinformatics/test.results/results/IFNy_A.fumigatusconidia_PBMC_7days/fold.1/IFNy_A.fumigatusconidia_PBMC_7days.rda")


x <- as.vector(as.matrix(IFNy_B.burgdorferi_PBMC_7days.coef[,c("coef.fold.1", "coef.fold.2", "coef.fold.3", "coef.fold.4", "coef.fold.5", "coef.fold.6","coef.fold.6", "coef.fold.7", "coef.fold.8", "coef.fold.9", "coef.fold.10")]))


ggplot2.violinplot(data=x)


p <- ggplot(data = NULL, aes(y = x, x="IFNy_B.burgdorferi_PBMC_7days")) +
  geom_violin(na.rm = T, fill="#1fa2ba") + 
  geom_boxplot(width = 0.2, fill="#42d9f4") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank())

p

