# libraries
library(data.table)
library(optparse)
library(glmnet)
library(MASS)
library(methods)
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/new.script/create.k.folds.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/new.script/cv.glmnet.wrapper.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/new.script/glmnet.wrapper.R")

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
phenotype <- opt$p
database <- opt$d

print("DEFINING VARIABLES")

path.to.cor <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/cor/cor.matrix/"
path.to.results <-"/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/10.imputation.10.fold.cor.threshold.cqtl.mapping."
molecular <- "/groups/umcg-wijmenga/tmp04/data/500FG/molecular_phenotypes/"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/"
cytokines <- "/groups/umcg-wijmenga/tmp04/data/500FG/cytokine_levels/500FG_log2_cytokine_levels.txt"
genetics <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/qtl/genetics.txt"
cqtl <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/qtl/"

# globals for testing
# cytokines <- "~/Bioinformatics/data/cytokines/cytokine.levels.test.txt"
# molecular <- "~/Bioinformatics/data/molecular.phenotypes/"
# gene.expression <- "~/Bioinformatics/data/cytokines/"
# path.to.cor <- "~/Bioinformatics/data/cor.matrix/"
# folds <- 1
# phenotype <- "TNFA_S.aureus_WB_48h"
# database <- "cytokines"
# path.to.results <- "~/Bioinformatics/results/"
# genetics <- "~/Bioinformatics/data/cqtl_mapping/500FG_dosage_matrix_cytokine_samples_only.txt"
# cqtl <- "~/Bioinformatics/data/cqtl_mapping/filtered_cQTL_pval_5e-4.txt"

# data
print("READING CELL COUNTS")
cell.counts <- read.table(paste0(molecular, "500FG_inverse_rank_normalized_cellcounts.txt"))
print("READING IMMUNOGLOBULIN LEVELS")
immunoglobulin.levels <- read.table(paste0(molecular, "500FG_log2_immunoglobulin_levels.txt"))
print("READING CYTOKINE LEVELS")
cytokines <- fread(cytokines); rownames(cytokines) <- cytokines$id; cytokines$id <- NULL

print("READING GENE EXPRESSION")
genes   <- fread(paste0(gene.expression, "500FG_RNAseq_TMM_normalized_read_counts.txt"), data.table = F); rownames(genes) <- genes$id; genes$id <- NULL
print("READING GENOTYPE DATA")
genetics <- fread(genetics, data.table = F); rownames(genetics) <- genetics$id; genetics$id <- NULL; genetics <- t(genetics)
print("READING QTL MAPPING")
cqtl <- fread(paste0(cqtl, database, ".txt"))


# for threshold
for(thr in c("5e-4", "5e-5", "5e-6")){
  # adjust path
  
  p <- paste0(path.to.results, thr, "/")
  
  # create dirs
  if(!dir.exists(file.path(p, database))){
    dir.create(file.path(p, database))
  }
  
  # create subdir
  if(!dir.exists(file.path(paste0(p, database, "/"), phenotype))){
    dir.create(file.path(paste0(p, database, "/"), phenotype))
  }
  

}

# assign database given
data <- get(database)
print("READ CORRELATION MATRIX")
# load cor
cor   <- fread(paste0(path.to.cor, database, ".cor.txt")); rownames(cor) <- colnames(data)

print("FILTERING GENES ON PHENOTYPE")
# filter genes per cytokine
tmp <- cor[which(rownames(cor) == phenotype), ]
g <- names(which(tmp[ , apply(tmp, 2, function(x) any(x < .05))] == TRUE))

# y
y <- data[[phenotype]]; names(y) <- rownames(data)

# x
print("REDUCING GENES ON SIGNIFICANCE")
genes <- genes[, g]

# reduce genetics on samples from phenotype
print("FILTERING GENOTYPE DATA ON PHENOTYPE")
int <- Reduce(intersect, list(names(y), rownames(genetics)))
genetics <- genetics[rownames(genetics) %in% int, ]




# threshold cgt
print("FOR THREE THRESHOLDS PERFORM ELASTIC NET")
for(thr in c("5e-6")){
  print(thr)
  print("FILTER QTL")
  cqtl.new <- cqtl[which(cqtl$pvalue < as.numeric(thr)), ]
  
  print("FILTER GENETICS")
  # intercept cqtl_mapping with genetics
  int <- Reduce(intersect, list(cqtl.new$snps, colnames(genetics)))
  genetics.new <- genetics[, colnames(genetics) %in% int ]
  
  print("FILTER GENES AND GENETICS")
  # intercept genetics with gene expression
  int <- Reduce(intersect, list(rownames(genetics.new), rownames(genes)))
  
  # bind frames together
  genes.new    <- genes[rownames(genes) %in% int, ]
  genetics.new <- genetics.new[rownames(genetics.new) %in% int, ]
  
  # new x
  print("COMBINE GENETICS AND EXPRESSION")
  x <- merge(x = genes.new, y = genetics.new, by = "row.names"); rownames(x) <- x$Row.names; x$Row.names <- NULL
  
  # print size
  print(length(colnames(x)))
  
  # Remove NA values
  print("REMOVING NA")
  y             <- y[!is.na(y)]
  x             <- na.omit(x)
  print("INTERSECT")
  ol            <- intersect(rownames(x), names(y))
  y             <- y[ol]
  x             <- x[ol, , drop=F]
  
  # perform elastic net
  for(folds in 1:10){
  
  # create dir for imputation if not exist
  if(!dir.exists(file.path(paste0(p, database, "/", phenotype, "/"), paste0("fold.", folds)))){
    dir.create(file.path(paste0(p, database, "/", phenotype, "/"), paste0("fold.", folds)))
  }
    
  print("PERFORM ELASTIC NET")
  result <- cv.glmnet.wrapper(y=y, x=x, kfold=10, remove.outliers=F)
  print(paste0("PERFORMED", folds))
  # store result
  save(result, file = paste0(path.to.results, thr, "/",
                             database, "/",
                             phenotype,
                             "/", "fold.", folds,
                             "/",
                             phenotype,
                             ".rda"))
  }
}
