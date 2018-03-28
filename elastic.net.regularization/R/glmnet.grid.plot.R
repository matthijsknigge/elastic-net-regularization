#' plot elastic net
#' @author Matthijs Knigge
#'
#' @param plot.list list, ggplot objects
#' @param x.lab vector character names of x-axis
#' @param title title of plot character
#' @param save boolean, save plot. default FASLE
#' @param return boolean, return instead of save. Default FALSE
#' @param path character to path to save, must include name of plot.
#' @param width numeric width of plot in cm
#' @param height numeric height of plot in cm
#' @keywords plot
#' @export
#' @examples
#' glmnet.grid.plot()
#' @return ggplot object
#'
#'
glmnet.grid.plot <- function(plot.list, x.lab, title, save = FALSE, return = FALSE, path, width = 10, height = 10){
  # libraries
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  require(GGally)
  
  # create grid
  pm <- ggmatrix(
    plot.list,
    nrow = 1, ncol = length(plot.list),
    xAxisLabels = x.lab,
    title = title
  )
  # save
  if(save){
    ggsave(plot = pm, filename = path, width = width, height = height, dpi = 300)
  }
  # return
  if(return){
    return(list(pm = pm))
  }
  
}
