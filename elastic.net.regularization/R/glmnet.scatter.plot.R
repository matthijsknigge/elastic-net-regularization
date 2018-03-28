#' Scatter plots
#' @author Matthijs Knigge
#'
#' @param x numeric vector. Mean correlation of interaction between immune phenotypes and genes
#' @param y numeric vector. Mean correlation of interaction between immune phenotypes and genes + genetics
#' @param title character. Title of the plot
#' @keywords linear regression
#' @export
#' @examples
#' glmnet.scatter.plot()
#' @return 
#' 
glmnet.scatter.plot <- function(x, y, title){
  # libraries
  require(ggplot2)
  
  # theme for plotting
  theme_set(theme_gray()) 
  
  # intercept x with y
  int <- Reduce(intersect, list(rownames(x), rownames(y)))
  # apply on x
  x <- x[rownames(x) %in% int,]
  # apply on y
  y <- y[rownames(y) %in% int,]
  
  # define frame
  tmp <- cbind(x, y); tmp <- as.data.frame(tmp)
  # remove NA
  tmp <- na.omit(tmp)
  
  # perform pairwise t.test
  tmp$z <- apply(tmp, 1, function(x) {
    tmp.x <- x[1:10]
    tmp.y <- x[11:20]
    
    t.test(tmp.x, tmp.y, paired = TRUE, alternative = "two.sided")$p.value
  })
  # mean X
  tmp$x.mean <- apply(tmp, 1, function(x) {
    tmp.x <- x[1:10]
    median(tmp.x)
  })
  # mean Y
  tmp$y.mean <- apply(tmp, 1, function(x) {
    tmp.y <- x[11:20]
    median(tmp.y)
  })
  # add color
  tmp$col <- ""
  # add color for x
  tmp[which(tmp$x.mean > tmp$y.mean), ]$col <- "purple"
  # add color for y
  tmp[which(tmp$y.mean > tmp$x.mean), ]$col <- "darkgreen"
  # add black
  tmp[which(tmp$z > 0.05), ]$col <- "black"
  
  # log scale
  tmp$z <- log10(tmp$z)*-1
  
  # initiate ggplot
  p <- ggplot(data = tmp, aes(x.mean, y.mean, size = z))
  # draw points
  p <- p + geom_point(color=tmp$col)
  # 45 angle line
  p <- p + geom_abline(intercept = 0, slope = 1, color = "red")
  
  # change name
  p <- p + guides(size=guide_legend(title="10LOG(T.TEST)"))
  
  # add labs
  p <- p + labs(subtitle=title,
                x="GENE EXPRESSION",
                y="GENE EXPRESSION + GENETICS",
                title="GENE EXPRESSION ~ GENE EXPRESSION + GENETICS ",
                caption = "MATTHIJS KNIGGE")
  # adjust theme
  p <- p + theme(axis.text=element_text(size=15),
                 axis.title=element_text(size=10,face="bold"),
                 plot.title = element_text(face="bold", size=20),
                 legend.position = "bottom")
  # plot
  return(list(plot = p))
}