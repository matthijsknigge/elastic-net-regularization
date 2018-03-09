#' plot elastic net
#' @author Matthijs Knigge
#'
#' @param obs vector float, containing predicted values for Y
#' @param pred vector float, containing observed values for X
#' @param fold vector character, containing fold identifier
#' @param identifier character name of phenotype
#' @param imp numeric, impuation identifier
#' @keywords plot
#' @export
#' @examples
#' glmnet.plot()
#' @return ggplot object
#' 
#'      
glmnet.plot <- function(obs, pred, fold, identifier, imp){
  # libraries
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  
  # theme for plotting
  theme_set(theme_gray()) 

  # set min limit for observation
  obs.min <- min(obs)
  # set max limit for observation
  obs.max <- max(obs)
  # set min limit for predicted
  pred.min <- min(pred)
  # set max limit for predicted
  pred.max <- max(pred)
  
  # initiaite ggplot object
  p <- ggplot(data = NULL, aes(x = obs, y = pred, color = fold))
  # add points
  p <- p + geom_point()
  # add loess fit
  p <- p + geom_smooth(method="lm", se=F)
  # add labs
  p <- p + labs(subtitle=paste0("10-FOLD 10-VALIDATION, IMPUTATION ", imp),
                y="PREDICTED",
                x="MEASURED",
                title=paste0("PREDICTED ~ MEASURED ", identifier),
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
  # add rug
  p <- p + geom_rug(data=NULL, aes(x=obs, y=pred, color=fold))
  
  # push to global scope
  obs  <<- obs
  pred <<- pred
  fold <<- fold
  
  # create density x-axis
  xdens <- axis_canvas(p, axis = "x")
  xdens <- xdens + geom_density(data = NULL, aes(x = obs, fill=fold, color=fold), alpha = .5, size = 0.2)
  
  # create density y-axis
  ydens <- axis_canvas(p, axis = "y", coord_flip = T)
  ydens <- ydens + geom_density(data = NULL, aes(x = pred, fill = fold, color = fold), alpha=.5, size=.2)
  ydens <- ydens + coord_flip() 
  
  #insert
  p <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
  p <- insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")
  
  # remove from global environment
  rm(obs, envir = .GlobalEnv); rm(pred, envir = .GlobalEnv); rm(fold, envir = .GlobalEnv)
  
  
  # save plot
  return(list(p = p))
 
}
