# libraries
library(data.table)

# functions
do.cor <- function(data.x, data.y, name){
  # get intersect
  print("finding intersept"); int <- Reduce(intersect, list(rownames(data.x), rownames(data.y)))
  # apply intersect
  print("applying to data.x"); data.x <- data.x[rownames(data.x) %in% int,]
  # filter genes
  print("applying to data.y"); data.y <- data.y[rownames(data.y) %in% int,]
  # make data.matrix
  print("converting data.x"); data.x <- data.matrix(data.x)
  # make data.matrix 
  print("converting data.y"); data.y <- data.matrix(data.y)
  
  # do cor
  print("predicting correlation"); cor.pvalues <- data.frame(apply(data.y, 2, function(x) { apply(data.x, 2,   
                                                                function(y) { cor.test(x,y, method = "spearman")$p.value }) })) 
  # add rownames
  cor.pvalues$identifier <- rownames(cor.pvalues)
  
  # write.output
  print("writing output"); write.table(x = cor.pvalues, file = paste0(path.to.output, name, ".txt"), quote = T, row.names = F, col.names = T)
}

# globals
molecular <- "/groups/umcg-wijmenga/tmp04/data/500FG/molecular_phenotypes/"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/"
path.to.output <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/cor/cor.matrix/"

# globals for testing
# molecular <- "Bioinformatics/data/molecular.phenotypes/"
# gene.expression <- "Bioinformatics/data/cytokines/"
# path.to.output  <- "Desktop/"

# info
# cell.counts.info <- read.table(paste0(molecular, "500FG_cellcounts_info.txt"))
# hormone.info <- read.table(paste0(molecular, "500FG_hormone_info.txt"))

# molecular phenotypes
arcsine <- read.table(paste0(molecular, "500FG_arcsine_transformed_taxonomy_abundance.txt")); rownames(arcsine) <- arcsine$X; arcsine$X <- NULL 
cell.counts<- read.table(paste0(molecular, "500FG_inverse_rank_normalized_cellcounts.txt"))
metabolites <- read.table(paste0(molecular, "500FG_log2_brainshake_metabolites.txt"))
protein.levels <- read.table(paste0(molecular, "500FG_log2_circulating_protein_levels.txt"))
hormone.levels <- read.table(paste0(molecular, "500FG_log2_hormone_levels.txt"))
immunoglobulin.levels <- read.table(paste0(molecular, "500FG_log2_immunoglobulin_levels.txt"))
platelet.counts <- read.table(paste0(molecular, "500FG_log2_platelet_count.txt"))
microbiome <- read.table(paste0(molecular, "500FG_microbiome_pathways.txt"))

# gene expression
genes   <- fread(paste0(gene.expression, "500FG_RNAseq_TMM_normalized_read_counts.txt")); rownames(genes) <- genes$id; genes$id <- NULL

# dor cor
print("arcsine.cor")
do.cor(data.x = arcsine, data.y = genes, name = "arcsine.cor")
print("cell.counts.cor")
do.cor(data.x = cell.counts, data.y = genes, name = "cell.counts.cor")
print("metabolites.cor")
do.cor(data.x = metabolites, data.y = genes, name = "metabolites.cor")
print("protein.levels.cor")
do.cor(data.x = protein.levels, data.y = genes, name = "protein.levels.cor")
print("hormone.levels.cor")
do.cor(data.x = hormone.levels, data.y = genes, name = "hormone.levels.cor")
print("immunoglobulin.levels.cor")
do.cor(data.x = immunoglobulin.levels, data.y = genes, name = "immunoglobulin.levels.cor")
print("platelet.counts.cor")
do.cor(data.x = platelet.counts, data.y = genes, name = "platelet.counts.cor")
print("microbiome.cor")
do.cor(data.x = microbiome, data.y = genes, name = "microbiome.cor")




