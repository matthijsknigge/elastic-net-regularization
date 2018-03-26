#' box plot elastic net
#' @author Matthijs Knigge
#'
#' @param coef vector float, correlation coefficients
#' @param group vector character of grouping
#' @param identifier vector character of phenotypes
#' @param title character title of plot
#' @param lab boolean, show mean correlation for boxplots
#' @keywords plot
#' @export
#' @examples
#' glmnet.box.plot()
#' @return ggplot object
#' 
#'      
glmnet.box.plot <- function(coef, group, identifier, title){
  # libraries
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  
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
  
  # create ggplot object
  p <- ggplot(data = tmp.df, aes(y = coef, x = identifier, fill = group, color = group))
  # create facet for each group
  p <- p + facet_grid(group~., scales = "free", space = "free")
  # add jitter points
  p <- p + geom_jitter(col="grey", alpha=.7, show.legend = FALSE)
  # add boxplot
  p <- p + geom_boxplot(width = 0.7, position=position_dodge(1), alpha=.7)
  # mean lines
  p <- p + geom_hline(data = df.mean.coef, aes(yintercept = mean, color = group), size=1, linetype="solid")
  # add zero line
  p <- p + geom_hline(yintercept = 0, size=.3, linetype="solid", color="black")
  # adjust theme
  p <- p + theme(
    axis.ticks.x=element_blank(), 
    axis.line.x = element_line(colour = 'white', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'white', size=0.5, linetype='solid'),
    legend.position = "bottom", 
    legend.title=element_blank())
  # add labels
  p <- p + labs(subtitle=title,
                x="IMMUNE PHENOTYPES",
                y="SPEARMEN COEFFICIENT",
                title="OVERVIEW PREDICTION POWER",
                caption = "MATTHIJS KNIGGE")
  # flip 
  p <- p + coord_flip()
  
  # save plot
  return(list(p = p))
  
}
