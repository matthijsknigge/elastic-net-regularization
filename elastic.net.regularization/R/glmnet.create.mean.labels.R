#' create mean labels
#' @author Matthijs Knigge
#'
#' @param data data.frame containing columns coef, group and identifier
#' @keywords labels
#' @export
#' @examples
#' glmnet.create.mean.labels()
#' @return ggplot object
#' 
#'      
glmnet.create.mean.labels <- function(data){
  # create new label
  for(id in unique(data$identifier)){
    # subset data
    tmp <- data[which(data$identifier == id), ]
    # calc mean of identifier
    m <- mean(tmp$coef, na.rm = T)
    # create new label
    lab <- paste0(id, "\t\t{", signif(x = m, digits = 3), "}")
    # insert into frame
    data[which(data$identifier == id), ]$identifier <- lab
  }
  return(list(d = data))
}