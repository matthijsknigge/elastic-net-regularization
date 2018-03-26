#' create mean labels
#' @author Matthijs Knigge
#'
#' @param identifier no comment yet
#' @param group no comment yet
#' @param coef no comment yet
#' @keywords labels
#' @export
#' @examples
#' glmnet.create.mean.labels()
#' @return ggplot object
#' 
#'      
glmnet.create.mean.labels <- function(identifier, group, coef){
  # create data.frame
  data <- data.frame(identifier = identifier, group = group, coef = coef, label = "", stringsAsFactors = F)
  # create new label
  for(id in unique(data$identifier)){
    # subset data
    tmp <- data[which(data$identifier == id), ]
    # calc mean of identifier
    m <- mean(tmp$coef, na.rm = T)
    # create new label
    lab <- paste0(id, " {", signif(x = m, digits = 3), "}")
    # insert into frame
    data[which(data$identifier == id), ]$label <- lab
  }
  # swap
  data$identifier <- data$label
  # remove
  data$label <- NULL
  # return
  return(list(d = data))
}