#' box plot elastic net
#' @author Olivier Bakker, Matthijs Knigge
#'
#' @param output.file character name of the outputfile
#' @param raw.cell.counts character absolute path to raw read counts
#' @keywords plot
#' @export
#' @examples
#' glmnet.TMM.normalize()
#' 
glmnet.TMM.normalize <- function(raw.cell.counts) {
  # libraries
  require(edgeR)
  require(data.table)

  # pre parse table to apptropriate format
  table           <- raw.cell.counts
  D               <- DGEList(counts=table)
  # Remove low expressed genes 500FG SPECIFIC

  # Remove zero genes
  keep          <- rowSums(D$counts == 0) < 0.8 * ncol(D$counts)
  D             <- D[keep, , keep.lib.sizes=FALSE]

  # Remove low varying genes
  var           <- apply(D$counts, 1, var)
  keep          <- var > 0.1
  D             <- D[keep, , keep.lib.sizes=FALSE]

  # TMM normilization
  d             <- calcNormFactors(D)
  scalar        <- d$samples$lib.size*d$samples$norm.factors/exp(mean(log(d$samples$lib.size*d$samples$norm.factors)))
  scal.mat      <- outer(rep(1,nrow(d$counts)), scalar)
  scaled.counts <- d$counts/scal.mat

  # return scaled counts
  return(list(scaled.counts = scaled.counts))
}







