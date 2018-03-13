#' plot elastic net
#' @author Matthijs Knigge
#'
#' @param group vector character of grouping
#' @param identifier vector character of phenotypes
#' @param title character title of plot
#' @param gene.amount.before numeric vector containing average genes before per imputation
#' @param gene.amount.after  numeric vector containing average genes after per imputation
#' @keywords plot
#' @export
#' @examples
#' glmnet.bar.plot()
#' @return ggplot object
#'
#'
glmnet.bar.plot <- function(identifier, group, gene.amount.before, gene.amount.after, title){
  # libraries
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  # bind data before
  data.before <- data.frame(identifier = identifier, group = group, gene.amount = gene.amount.before, when = "before")
  # bind data after
  data.after <- data.frame(identifier = identifier, group = group, gene.amount = gene.amount.after, when = "after")
  # bind together
  data.tot <- rbind(data.before, data.after)

  # set theme
  theme_set(theme_gray())
  # initiate ggplot
  p <- ggplot(data = data.tot, aes(x = identifier, y=gene.amount, fill = when, color = when))
  # create grid
  p <- p + facet_grid(group~., scales = "free", space = "free")
  # bars
  p <- p + geom_bar(stat="identity", width = .4)
  # set theme
  p <- p + theme(
    axis.ticks.x=element_blank(),
    axis.line.x = element_line(colour = 'white', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'white', size=0.5, linetype='solid'),
    legend.position = "bottom",
    legend.title=element_blank())
  # add labels
  p <- p + labs(subtitle=title,
                x="IMMUNE PHENOTYPES",
                y="AMOUNT OF GENES",
                title="AVERAGE AMOUNT GENES IN PREDICTED MODELS",
                caption = "MATTHIJS KNIGGE")
  # flip y and x
  p <- p + coord_flip()


  # save plot
  return(list(p = p))

}
