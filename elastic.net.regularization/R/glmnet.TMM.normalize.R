#' TMM normalize read counts
#'
#' @return 
#' @export
#' @examples
#' glmnet.TMM.normalize()
glmnet.TMM.normalize <- function() {
  args <- commandArgs(TRUE)
  
  require(edgeR)
  require(data.table)
  install.packages("edgeR")
  
  countsFile=args[1]
  normOut=args[2]
  
  table           <- fread(countsFile, data.table=F)
  rownames(table) <- table[,1]
  table           <- table[,-1]
  D               <- DGEList(counts=table)
  
  cat("[INFO] ", nrow(D), " genes in input\n")
  cat("[INFO] ", ncol(D), " samples in input\n")
  
  # Remove low expressed genes 500FG SPECIFIC
  
  # Remove zero genes
  keep          <- rowSums(D$counts == 0) < 0.8 * ncol(D$counts)
  D             <- D[keep, , keep.lib.sizes=FALSE]
  cat("[INFO] Removed genes with more then 80% zeroes.\n")
  
  # Remove low varying genes
  var           <- apply(D$counts, 1, var)
  keep          <- var > 0.1
  D             <- D[keep, , keep.lib.sizes=FALSE]
  
  cat("[INFO] ", nrow(D), " genes remaining after filtering.\n")
  
  
  # TMM normilization
  d             <- calcNormFactors(D)
  scalar        <- d$samples$lib.size*d$samples$norm.factors/exp(mean(log(d$samples$lib.size*d$samples$norm.factors)))
  scal.mat      <- outer(rep(1,nrow(d$counts)), scalar)
  scaled.counts <- d$counts/scal.mat
  
  write.table(t(scaled.counts), file = normOut, sep = "\t", row.names=TRUE, quote=FALSE, col.names = NA)
  cat("[INFO] Output written.\n")
  
}














fwrite(col.names = T, row.names = F)