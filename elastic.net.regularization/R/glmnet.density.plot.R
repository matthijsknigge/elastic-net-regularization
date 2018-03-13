#' density plot elastic net
#' @author Matthijs Knigge
#'
#' @param coef vector float, correlation coefficients
#' @param group vector character of grouping
#' @param identifier vector character of phenotypes
#' @param title character title of plot
#' @keywords plot
#' @export
#' @examples
#' glmnet.density.plot()
#' @return ggplot object
#'
#'
glmnet.density.plot <- function(coef, group, identifier, title){
  # libraries
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  require(ggridges)

  # set theme
  theme_set(theme_gray())

  # create tmp.df for calculating mean
  tmp.df <- data.frame(coef = coef, group = group, identifier = identifier)
  # create df.mean.coef for storing mean
  df.mean.coef <- data.frame(identifier = character(), mean = numeric())
  # fill df.mean.coef
  for(p in unique(tmp.df$group)){
    # split frame for unique identifiers
    tmp <- tmp.df[which(tmp.df$group == p), ]
    # identifier
    id <- p
    # mean
    mean <- mean(tmp$coef, na.rm = T)
    # bind to df.mean.coef
    df.mean.coef <- rbind(df.mean.coef, data.frame(group = id,
                                                   mean = mean))
  }
  # add levels to df.mean.coef
  df.mean.coef$group <- factor(df.mean.coef$group)
  # add levels to tmp.df
  tmp.df$group <- factor(tmp.df$group)


  # initiate ggplot object
  p <- ggplot(data = tmp.df, aes(y = identifier, x = coef, fill = group, color = group))
  # facets
  p <- p + facet_grid(group~., scales = "free", space = "free")
  # points
  p <- p + geom_jitter(col="grey", alpha=.7, show.legend = FALSE)
  # density
  p <- p + geom_density_ridges(alpha = 0.6, scale = .6)
  # vertical line of mean power
  p <- p + geom_vline(data = df.mean.coef, aes(xintercept = mean, color = group), size=1, linetype="solid")
  # base pine
  p <- p + geom_vline(xintercept = 0, color = "black", size=.3, linetype="solid")
  # set theme
  p <- p + theme(
    axis.ticks.x=element_blank(),
    axis.line.x = element_line(colour = 'white', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'white', size=0.5, linetype='solid'),
    legend.position = "bottom",
    legend.title=element_blank())
  # add labs
  p <- p + labs(subtitle=title,
                    y="IMMUNE PHENOTYPES",
                    x="DENSITY",
                    title="DENSITY PREDICTION POWER",
                    caption = "MATTHIJS KNIGGE")
  
  # save plot
  return(list(p = p))

}
