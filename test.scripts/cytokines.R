# libraries
library(ggplot2)
library(optparse)
library(data.table)
library(mgcv)
library(Matrix)
library(MASS)  
library(glmnet)
library(methods)

# options
option_list = list(
  make_option(c("-c", "--cytokine"), 
              help=""),
  make_option(c("-f", "--fold"), 
              help="")
)
opt_parser    <- OptionParser(option_list=option_list, description="")
opt           <- parse_args(opt_parser)

# globals
folds <- opt$f
cytokine.name <- opt$c
cytokine.levels <- "/groups/umcg-wijmenga/tmp04/data/500FG/cytokine_levels/500FG_log2_cytokine_levels.txt"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/500FG_RNAseq_TMM_normalized_read_counts.txt"
path.to.cor <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/cor/cor.matrix.txt"
path.to.results <-"/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/"
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/create.k.folds.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/cv.glmnet.wrapper.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/glmnet.wrapper.R")

# create dir for cytokine if not exist
if(!dir.exists(file.path("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/", cytokine.name))){
  dir.create(file.path("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/", cytokine.name))
}

# create dir for imputation if not exist
if(!dir.exists(file.path(paste0("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/", cytokine.name, "/"), paste0("fold.", folds)))){
  dir.create(file.path(paste0("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/", cytokine.name, "/"), paste0("fold.", folds)))
}

# globals for testing
# cytokine.levels <- "~/Bioinformatics/data/cytokine.levels.test.txt"
# folds <- 10
# cytokine.levels.info <- "~/Bioinformatics/data/cytokine.levels.info.txt"
# gene.expression <- "~/Bioinformatics/data/expression.test.txt"
# cytokine.name <- "TNFA_S.typhimurium_macroPG_24h"
# path.to.results <- "~/Bioinformatics/test.results/"
# source("Bioinformatics/interplay.gut.microbiome.immune.system/test.scripts/create.k.folds.R")
# source("Bioinformatics/interplay.gut.microbiome.immune.system/test.scripts/cv.glmnet.wrapper.R")
# source("Bioinformatics/interplay.gut.microbiome.immune.system/test.scripts/glmnet.wrapper.R")
# path.to.results <- "~/Bioinformatics/test.results/"
# path.to.cor <- "~/Bioinformatics/data/cor.matrix.txt" 

# data
c.l   <- fread(cytokine.levels); rownames(c.l) <- c.l$id; c.l$id <- NULL
g.e   <- fread(gene.expression, data.table = T); rownames(g.e) <- g.e$id; g.e$id <- NULL; g.e <- data.matrix(g.e)
cor   <- fread(path.to.cor); rownames(cor) <- cor$cytokines; cor$cytokines <- NULL

# filter genes per cytokine
tmp <- cor[which(rownames(cor) == cytokine.name), ]
genes <- names(which(tmp[ , apply(tmp, 2, function(x) any(x < .05))] == TRUE))

# y
y <- c.l[[cytokine.name]]; names(y) <- rownames(c.l)
# x
x <- g.e[, genes]

result <- cv.glmnet.wrapper(y=y, x=x, kfold=10, remove.outliers=F)
# store result
save(result, file = paste0(path.to.results, 
                           cytokine.name,
                           "/", "fold.", folds,
                           "/", 
                           cytokine.name,
                           ".rda"))




