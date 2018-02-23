# libraries
library(ggplot2)
library(cross.validation)
library(elastic.net.regularization)
library(optparse)
library(data.table)
library(mgcv)


# options
option_list = list(
  make_option(c("-c", "--cytokine"),
              help=" "))

opt_parser    <- OptionParser(option_list=option_list, description="")
opt           <- parse_args(opt_parser)


# globals
cytokine.name <- opt$c
cytokine.levels <- "/groups/umcg-wijmenga/tmp04/data/500FG/cytokine_levels/500FG_log2_cytokine_levels.txt"
cytokine.levels.info <- "/groups/umcg-wijmenga/tmp04/data/500FG/cytokine_levels/500FG_cytokine_info.txt"
gene.expression <- "/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/500FG_RNAseq_TMM_normalized_read_counts.txt"
path.to.results <-"/groups/umcg-wijmenga/tmp04/elastic.net/results/"

# globals for testing
# cytokine.levels <- "~/Bioinformatics/data/cytokine.levels.test.txt"
# cytokine.levels.info <- "~/Bioinformatics/data/cytokine.levels.info.txt"
# gene.expression <- "~/Bioinformatics/data/expression.test.txt"
# cytokine.name <- "IFNy_LPS_WB_48h"

# data
c.l.i <- fread(cytokine.levels.info)
c.l   <- fread(cytokine.levels); rownames(c.l) <- c.l$id; c.l$id <- NULL
g.e   <- fread(gene.expression, stringsAsFactors = F, data.table = T); rownames(g.e) <- g.e$id; g.e$id <- NULL; g.e <- data.matrix(g.e)

# y
y <- c.l[[cytokine.name]]; names(y) <- rownames(c.l)
# x
x <- g.e

# start performing elastic net
result <- cv.glmnet.wrapper(y=y, x=x, kfold=10)

# store result
saveRDS(object = result, file = paste0(path.to.results, cytokine.name, ".rda"))
IFNy_LPS_WB_48h <- mod2 <- readRDS(paste0(path.to.results, cytokine.name, ".rda"))

devtools::install_bitbucket("matthijsknigge/Interplay.Between.Gut.Microbiome.and.Immune.System/cross.validation")
