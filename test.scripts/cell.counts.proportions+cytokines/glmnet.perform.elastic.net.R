# libraries
library(data.table)
library(MASS)
library(methods)
library(glmnet)
library(optparse)
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/cell.counts.proportion/create.k.folds.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/cell.counts.proportion/cv.glmnet.wrapper.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/cell.counts.proportion/glmnet.wrapper.R")


# options
option_list = list(
  make_option(c("-p", "--phenotype"),
              help=""),
  make_option(c("-d", "--database"),
              help="")
)
opt_parser    <- OptionParser(option_list=option_list, description="")
opt           <- parse_args(opt_parser)

# globals
folds <- 1
phenotype <- opt$p
database <- opt$d

cell.counts.proportion <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/data/phenotypes/cell.counts.proportion.txt"
path.to.cor <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/data/cor.matrix/"
path.to.results <-"/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/1.imputation.10.fold.cor.threshold.proportion.cell.counts/"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/gene.expression/"
molecular <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/data/phenotypes"


# globals for testing
# source("~/Bioinformatics/interplay.gut.microbiome.immune.system/cross.validation/R/create.k.folds.R")
# source("~/Bioinformatics/interplay.gut.microbiome.immune.system/elastic.net.regularization/R/glmnet.wrapper.R")
# source("~/Bioinformatics/interplay.gut.microbiome.immune.system/elastic.net.regularization/R/cv.glmnet.wrapper.R")
# phenotype <- "IT11"
# database <- "cell.counts.proportion"
# gene.expression <- "~/Bioinformatics/data/cytokines/"
# path.to.cor <- "~/Bioinformatics/data/cor.matrix/"
# cell.counts.proportion <- "~/Bioinformatics/data/molecular.phenotypes/500FG_proportion_cellcounts.txt"

# load data
cell.counts <- read.table(paste0(molecular, "500FG_inverse_rank_normalized_cellcounts.txt"))
immunoglobulin.levels <- read.table(paste0(molecular, "500FG_log2_immunoglobulin_levels.txt"))
cytokines <- fread(cytokines); rownames(cytokines) <- cytokines$id; cytokines$id <- NULL
cell.counts.proportion <- read.table(cell.counts.proportion)
genes   <- fread(paste0(gene.expression, "500FG_RNAseq_TMM_normalized_read_counts.txt")); rownames(genes) <- genes$id; genes$id <- NULL; genes <- data.matrix(genes)

# create dir
if(!dir.exists(file.path(paste0(path.to.results), phenotype))){
  dir.create(file.path(paste0(path.to.results), phenotype))
}

# assign database given
data <- get(database)

# load cor
cor   <- fread(paste0(path.to.cor, database, ".cor.txt")); rownames(cor) <- colnames(data)

# y
y <- data[[phenotype]]; names(y) <- rownames(data)

# filter genes per cytokine
tmp <- cor[which(rownames(cor) == phenotype), ]
g <- names(which(tmp[ , apply(tmp, 2, function(x) any(x < .05))] == TRUE))


# x
x <- genes[, g]

# for folds
for(i in 1:folds){
  
  # create dir
  dir.create(file.path(paste0(path.to.results, phenotype, "/", "fold.", i)))

  # perform elastic net
  result <- cv.glmnet.wrapper(y=y, x=x, kfold=10, remove.outliers=F)
  
  # store result
  save(result, file = paste0(path.to.results,
                             phenotype,
                             "/", "fold.", i,
                             "/",
                             phenotype,
                             ".rda"))
}






