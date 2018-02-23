# libraries
library(data.table)
library(optparse)
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/create.k.folds.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/cv.glmnet.wrapper.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/glmnet.wrapper.R")

# options
option_list = list(
  make_option(c("-p", "--phenotype"), 
              help=""),
  make_option(c("-f", "--fold"), 
              help=""),
  make_option(c("-d", "--database"), 
              help="")
)
opt_parser    <- OptionParser(option_list=option_list, description="")
opt           <- parse_args(opt_parser)

# globals
folds <- opt$f
phenotype <- opt$p
database <- opt$d

path.to.cor <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/cor/cor.matrix/"
path.to.results <-"/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/"
molecular <- "/groups/umcg-wijmenga/tmp04/data/500FG/molecular_phenotypes/"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/"
cytokines <- "/groups/umcg-wijmenga/tmp04/data/500FG/cytokine_levels/500FG_log2_cytokine_levels.txt"

# globals for testing
# cytokines <- "Bioinformatics/data/cytokines/cytokine.levels.test.txt"
# molecular <- "Bioinformatics/data/molecular.phenotypes/"
# gene.expression <- "Bioinformatics/data/cytokines/"
# path.to.cor <- "Bioinformatics/data/cor.matrix/"
# folds <- 1
# phenotype <- "PWY.5508..adenosylcobalamin.biosynthesis.from.cobyrinate.a.c.diamide.II"
# database <- "microbiome"
# path.to.results <- "Bioinformatics/results/"

# data
arcsine <- read.table(paste0(molecular, "500FG_arcsine_transformed_taxonomy_abundance.txt")); rownames(arcsine) <- arcsine$X; arcsine$X <- NULL 
cell.counts<- read.table(paste0(molecular, "500FG_inverse_rank_normalized_cellcounts.txt"))
metabolites <- read.table(paste0(molecular, "500FG_log2_brainshake_metabolites.txt"))
protein.levels <- read.table(paste0(molecular, "500FG_log2_circulating_protein_levels.txt"))
hormone.levels <- read.table(paste0(molecular, "500FG_log2_hormone_levels.txt"))
immunoglobulin.levels <- read.table(paste0(molecular, "500FG_log2_immunoglobulin_levels.txt"))
platelet.counts <- read.table(paste0(molecular, "500FG_log2_platelet_count.txt"))
microbiome <- read.table(paste0(molecular, "500FG_microbiome_pathways.txt"))
cytokines <- fread(cytokines); rownames(cytokines) <- cytokines$id; cytokines$id <- NULL
genes   <- fread(paste0(gene.expression, "500FG_RNAseq_TMM_normalized_read_counts.txt")); rownames(genes) <- genes$id; genes$id <- NULL


# create dirs
if(!dir.exists(file.path(path.to.results, database))){
  dir.create(file.path(path.to.results, database))
}

# create subdir
if(!dir.exists(file.path(paste0(path.to.results, database, "/"), phenotype))){
  dir.create(file.path(paste0(path.to.results, database, "/"), phenotype))
}

# create dir for imputation if not exist
if(!dir.exists(file.path(paste0(path.to.results, database, "/", phenotype, "/"), paste0("fold.", folds)))){
  dir.create(file.path(paste0(path.to.results, database, "/", phenotype, "/"), paste0("fold.", folds)))
}

# load cor
cor   <- fread(paste0(path.to.cor, database, ".cor.txt")); rownames(cor) <- cor$identifier; cor$identifier <- NULL

# filter genes per cytokine
tmp <- cor[which(rownames(cor) == phenotype), ]
g <- names(which(tmp[ , apply(tmp, 2, function(x) any(x < .05))] == TRUE))

# assign database given 
data <- get(database)

# y
y <- data[[phenotype]]; names(y) <- rownames(data)

# x
x <- genes[, g]

# perform elastic net
result <- cv.glmnet.wrapper(y=y, x=x, kfold=10, remove.outliers=F)

# store result
save(result, file = paste0(path.to.results, 
                           database,
                           phenotype,
                           "/", "fold.", folds,
                           "/", 
                           phenotype,
                           ".rda"))
# cytokine IFNy_LPS_WB_48h

