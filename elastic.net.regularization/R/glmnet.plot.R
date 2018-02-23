#' plot elastic net
#' @author Matthijs Knigge
#'
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @keywords plot
#' @export
#' @examples
#' glmnet.plot()
#' @return ggplot object
#' 
#'      
glmnet.plot <- function(glmnet.data.frame, ylab, xlab, title, interval, filename){
  # libraries
  require(ggplot2)
  require(gridExtra)

    # initiaite ggplot object
  p <- ggplot(data = imp.sub, aes(x = obs, y = pred, color = fold))
  # add points
  p <- p + geom_point()
  # add loess fit
  p <- p + geom_smooth(method="lm", se=F)
  # add labs
  p <- p + labs(subtitle="10-FOLD 10-VALIDATION",
                y="PREDICTED",
                x="MEASURED",
                title=paste0("PREDICTED ~ MEASURED ", c,"\n", imp),
                caption = "MATTHIJS KNIGGE")
  # adjust theme
  p <- p + theme(axis.text=element_text(size=15),
                 axis.title=element_text(size=10,face="bold"),
                 plot.title = element_text(face="bold", size=20),
                 legend.position = "bottom")
  # rescale x-axis for animation
  p <- p + scale_x_continuous(limits=c(obs.min, obs.max))
  # rescale y-axis for animation
  p <- p + scale_y_continuous(limits=c(pred.min, pred.max))
  
  # create density x-axis
  xdens <- axis_canvas(p, axis = "x")
  xdens <- xdens + geom_density(data = imp.sub, aes(x = obs, fill=fold, color=fold), alpha = .5, size = 0.2)
  
  # create density y-axis
  ydens <- axis_canvas(p, axis = "y", coord_flip = T)
  ydens <- ydens + geom_density(data = imp.sub, aes(x = pred, fill = fold, color = fold), alpha=.5, size=.2)
  ydens <- ydens + coord_flip() 
  
  #insert
  p <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  p <- insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")
  
  # save plot
  ggsave(filename = paste0(path.to.output, c, "/", c, ".", imp, ".png"), plot = p,
         width = 10, height = 10, dpi = 150)
 
}
