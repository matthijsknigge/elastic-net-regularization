#' Scatter plots
#' @author Matthijs Knigge
#'
#' @param x numeric vector. Mean correlation of interaction between immune phenotypes and genes
#' @param y numeric vector. Mean correlation of interaction between immune phenotypes and genes + genetics
#' @param title character. Title of plot
#' @param sub.itle character. Title of the plot
#' @keywords linear regression
#' @export
#' @examples
#' glmnet.scatter.plot()
#' @return 
#' 
glmnet.scatter.plot <- function(x, y, title, sub.title){
  # libraries
  require(ggplot2)

  # theme for plotting
  theme_set(theme_gray()) 
  
  # remove NA's
  x <- na.omit(x)
  # remove NA's
  y <- na.omit(y)
  
  # intercept
  int <- Reduce(intersect, list(unique(x$identifier), unique(y$identifier)))
  
  # apply int on x
  x <- x[x$identifier %in% int, ]
  
  # apply int on y
  y <- y[y$identifier %in% int, ]
  
  # create new.frame
  tmp.data <<- data.frame(x.mean = numeric(), y.mean = numeric(), t.test = numeric())
  
  # for every phenotype
  for(id in int){
    # get unique x
    tmp.x <- x[which(x$identifier == id), ]
    # get unique y
    tmp.y <- y[which(y$identifier == id), ]
    # test for observations
    if(nrow(tmp.x) != 1 & nrow(tmp.y) != 1){
      # median x
      median.x <- median(tmp.x$coef)
      # median y
      median.y <- median(tmp.y$coef)
      # test
      tmp.t.test <- t.test(tmp.x$coef, tmp.y$coef, paired = FALSE, alternative = "two.sided")$p.value
      # bind to frame
      tmp.data <<- rbind(tmp.data, data.frame("median.x" = median.x, "median.y" = median.y, "t.test" = tmp.t.test))
    }
  }
  
  # add color
  tmp.data$col <- ""
  
  # if possible
  if(nrow(tmp.data[which(tmp.data$median.x > tmp.data$median.y), ]) > 0){
    # add color for x
    tmp.data[which(tmp.data$median.x > tmp.data$median.y), ]$col <- "orange"
  }

  # if possible
  if(nrow(tmp.data[which(tmp.data$median.y > tmp.data$median.x), ]) > 0){
    # add color for y
    tmp.data[which(tmp.data$median.y > tmp.data$median.x), ]$col <- "blue"
  }
  
  # add black
  if(nrow(tmp.data[which(tmp.data$t.test > 0.05), ]) > 0){
    tmp.data[which(tmp.data$t.test > 0.05), ]$col <- "black"
  }
  
  # log scale
  tmp.data$t.test <- log10(tmp.data$t.test)*-1
  
  # initiate ggplot
  p <- ggplot(data = tmp.data, aes(x = median.x, y = median.y, size = t.test))
  # draw points
  p <- p + geom_point(color=tmp.data$col)
  # 45 angle line
  p <- p + geom_abline(intercept = 0, slope = 1, color = "red")
  # add rug
  # add rug
  p <- p + geom_rug(data=NULL, x=tmp.data$median.x, y=tmp.data$median.y, color=tmp.data$col, size=.8)
  
  # change name
  p <- p + guides(size=guide_legend(title="10LOG(T.TEST)"))
  # xlim
  p <- p + xlim(range(tmp.data$median.x, tmp.data$median.y))
  # ylim
  p <- p + ylim(range(tmp.data$median.y, tmp.data$median.x))
  # add labs
  p <- p + labs(subtitle=sub.title,
                x="GENE EXPRESSION",
                y="GENE EXPRESSION + GENETICS",
                title=paste0("GENE EXPRESSION ~ GENE EXPRESSION + GENETICS ", "[", title, "]"),
                caption = "MATTHIJS KNIGGE")
  # adjust theme
  p <- p + theme(axis.text=element_text(size=15),
                 axis.title=element_text(size=10,face="bold"),
                 plot.title = element_text(face="bold", size=20),
                 legend.position = "bottom")
  
  
    # plot
  return(list(plot = p))
}