#' calculate correlation matrix
#' @author Matthijs Knigge
#'
#' @param data.x matrix n*m, must have row, and colnames.
#' @param data.y matrix n*m, must have row, and colnames.
#' @param write  boolean, default = FALSE. Set TRUE for storing matrix.
#' @param path   character, default = NULL. Path to save matrix.
#' @param melt   boolean, default = FALSE. Set TRUE, if need to return or write melted version of correlation.
#' @keywords cor correlation
#' @export
#' @examples
#' glmnet.cor()
#' @return correlation matrix between y and x 
#' 
#'         
glmnet.cor <- function(data.x, data.y, write = FALSE, path = NULL, melt = FALSE, method = c("spearman")){
  # libraries
  
  # get intersect
  int <- Reduce(intersect, list(rownames(data.x), rownames(data.y)))
  # apply intersect
  data.x <- data.x[rownames(data.x) %in% int,]
  # filter genes
  data.y <- data.y[rownames(data.y) %in% int,]
  # make data.matrix
  data.x <- data.matrix(data.x)
  # make data.matrix 
  data.y <- data.matrix(data.y)
  # do cor
  cor.matrix <- data.frame(apply(data.y, 2, function(x) { apply(data.x, 2,   
                          function(y) { cor.test(x,y, method = "spearman")$p.value }) })) 
  # add rownames
  cor.matrix$identifier <- rownames(cor.matrix)
  
  # write.output
  if(write){
    # if path is !NULL
    if(!is.null(path)){
      write.table(x = cor.matrix, file = paste0(path, ".txt"), quote = T, row.names = F, col.names = T)
    }
  }
  # return
  return(list(cor.matrix = cor.matrix))
}
1121*100
