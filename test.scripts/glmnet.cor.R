# libraries
library(data.table)

# globals

# cytokines
print("reading cytokines")
cytokines <- fread("/groups/umcg-wijmenga/tmp04/data/500FG/cytokine_levels/500FG_log2_cytokine_levels.txt")
rownames(cytokines) <- cytokines$id; cytokines$id <- NULL
print("done")

# genes
print("reading genes")
genes     <- fread("/groups/umcg-wijmenga/tmp04/data/500FG/gene_expression/500FG_RNAseq_TMM_normalized_read_counts.txt")
rownames(genes) <- genes$id; genes$id <- NULL
print("done")

# first intersect the data
print("find intercept between cytonines and genes")
int <- Reduce(intersect, list(rownames(cytokines), rownames(genes)))
print("apply intersect")
# filter cytokines
cytokines <- cytokines[rownames(cytokines) %in% int,]
# filter genes
genes     <- genes[rownames(genes) %in% int,]

# make data.matrix for genes
genes <- data.matrix(genes)
# make data.matrix for cytokines
cytokines <- data.matrix(cytokines)

print("start predicting correlations")
# functions
cor.pvalues <- data.frame(apply(genes, 2, function(x) { apply(cytokines, 2,   
                      function(y) { cor.test(x,y, method = "spearman")$p.value }) }))

cor.pvalues$cytokines <- rownames(cor.pvalues)

# write output
write.table(x = cor.pvalues, file = "/groups/umcg-wijmenga/tmp04/umcg-mknigge/elastic.net/cor/cor.matrix.txt", quote = T, row.names = F, col.names = T)



