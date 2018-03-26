#' plot elastic net
#' @author Matthijs Knigge
#'
#' @param p.1 ggplot object
#' @param p.2 ggplot object
#' @param p.3 ggplot object
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
glmnet.grid.plot <- function(p.1, p.2, p.3, p.4, p.5, x.lab, title, save = FALSE, return = FALSE, path, width = 10, height = 10){
  # libraries
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  require(GGally)
  
  # create grid
  pm <- ggmatrix(
    list(p.1, p.2, p.3, p.4, p.5),
    nrow = 1, ncol = 5,
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
