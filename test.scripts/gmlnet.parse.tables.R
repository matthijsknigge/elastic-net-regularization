# libraries
library(data.table)
library(optparse)
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/create.k.folds.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/cv.glmnet.wrapper.R")
source("/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/script/glmnet.wrapper.R")

# options
option_list = list(
  make_option(c("-c", "--cytokine"), 
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
cytokine.name <- opt$c
database <- opt$d

path.to.cor <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/cor/cor.matrix/"
path.to.results <-"/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/results/"
molecular <- "/groups/umcg-wijmenga/tmp04/data/500FG/molecular_phenotypes/"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/"


# globals for testing
# molecular <- "Bioinformatics/data/molecular.phenotypes/"
# gene.expression <- "Bioinformatics/data/cytokines/"

# data
arcsine <- read.table(paste0(molecular, "500FG_arcsine_transformed_taxonomy_abundance.txt")); rownames(arcsine) <- arcsine$X; arcsine$X <- NULL 
cell.counts<- read.table(paste0(molecular, "500FG_inverse_rank_normalized_cellcounts.txt"))
metabolites <- read.table(paste0(molecular, "500FG_log2_brainshake_metabolites.txt"))
protein.levels <- read.table(paste0(molecular, "500FG_log2_circulating_protein_levels.txt"))
hormone.levels <- read.table(paste0(molecular, "500FG_log2_hormone_levels.txt"))
immunoglobulin.levels <- read.table(paste0(molecular, "500FG_log2_immunoglobulin_levels.txt"))
platelet.counts <- read.table(paste0(molecular, "500FG_log2_platelet_count.txt"))
microbiome <- read.table(paste0(molecular, "500FG_microbiome_pathways.txt"))

# functions

# scripts



