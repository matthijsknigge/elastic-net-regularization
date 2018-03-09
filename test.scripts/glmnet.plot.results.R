# libraries
library(ggplot2)
library(data.table)
library(gridExtra)
library(cowplot)
library(data.table)
library("ggsci")
library("ggridges")

# global variables
path.to.output <- "~/Bioinformatics/test.results/plots/"

# data
obs.pred <- fread("~/Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/data.summary/df.pred.obs.txt", col.names = c("sample", "pred", "obs", "fold", "imputation", "cytokine", "stimulation"))
coef <- fread("~/Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/data.summary/df.coef.txt", col.names = c("imputation", "cytokine", "stimulation","coef.fold.1", "coef.fold.2", "coef.fold.3", "coef.fold.4","coef.fold.5", "coef.fold.6", "coef.fold.7", "coef.fold.8","coef.fold.9", "coef.fold.10", "amount.genes"))
cor <- fread("~/Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/data.summary/cor.matrix.txt"); 
coef.new <- fread("Bioinformatics/results/10.imputation.10.fold.immune.phenotypes/data.summary/df.coef.molten.txt")
cytokine.info <- fread("Bioinformatics/data/cytokine.levels.info.txt")
cor.molt <- fread("Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/data.summary/cor.molten.txt")
coef.molt <- fread("Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/data.summary/df.coef.molten.txt")
c.mean <- fread("Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/data.summary/cytokines.mean.genes.txt")


# functions
# theme
theme_set(theme_gray()) 

# Load big5 dataset:

# per cytokine
for(c in unique(obs.pred$cytokine)){
  # subset matrix on immune phenotype
  tmp <- obs.pred[which(obs.pred$cytokine == c), ]
  # set min limit for observation
  obs.min <- min(tmp$obs)
  # set max limit for observation
  obs.max <- max(tmp$obs)
  # set min limit for predicted
  pred.min <- min(tmp$pred)
  # set max limit for predicted
  pred.max <- max(tmp$pred)
  
  # per impuation
  for(imp in unique(tmp$imputation)){
    
    
    
    # subset on imputation
    imp.sub <- tmp[which(tmp$imputation == imp), ]
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
}




coef.new$stimulation <- factor(coef.new$stimulation)
unique(coef.new$stimulation)
dummy2 <- data.frame(stimulation = c("IFNy", "IL17", "IL1b", "IL22", "IL6", "TNFA"), sti = c(0.116673,-0.07609174,
                                                                                             0.2111578,0.1997144,
                                                                                             0.1221601,0.1905662))
dummy2$stimulation <- factor(dummy2$stimulation)




p <- ggplot(data = coef.new, aes(y = coef, x = cytokine, fill = stimulation, color = stimulation))
p <- p + facet_grid(stimulation~., scales = "free", space = "free")
p <- p + geom_boxplot(width = 0.7, position=position_dodge(1), alpha=.7)
p <- p + geom_jitter(col="grey", alpha=.7, show.legend = FALSE)
p <- p + geom_hline(data = dummy2, aes(yintercept = sti, color = stimulation), size=1, linetype="solid")
p <- p + geom_hline(yintercept = 0, size=.3, linetype="solid", color="black")
p <- p + theme(
               axis.ticks.x=element_blank(), 
               axis.line.x = element_line(colour = 'white', size=0.5, linetype='solid'),
               axis.line.y = element_line(colour = 'white', size=0.5, linetype='solid'),
               legend.position = "bottom", 
               legend.title=element_blank())
p <- p + labs(subtitle="10-FOLD 10-VALIDATION ON 91 IMMUNE PHENOTYPES",
             x="IMMUNE PHENOTYPES",
             y="SPEARMEN COEFFICIENT",
             title="OVERVIEW PREDICTION POWER",
             caption = "MATTHIJS KNIGGE")
p <- p + coord_flip()

p
ggsave(plot = p, filename = "Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/overview/prediction.power.png", height = 15, width = 10, dpi = 300)



p.1 <- ggplot(data = coef.new, aes(y = cytokine, x = coef, fill = stimulation, color = stimulation))

p.1 <- p.1 + facet_grid(stimulation~., scales = "free", space = "free")
p.1 <- p.1 + geom_jitter(col="grey", alpha=.7, show.legend = FALSE)
p.1 <- p.1 + geom_density_ridges(alpha = 0.6, scale = .6) 
p.1 <- p.1 + geom_vline(data = dummy2, aes(xintercept = sti, color = stimulation), size=1, linetype="solid")
p.1 <- p.1 + geom_vline(xintercept = 0, color = "black", size=.3, linetype="solid")
p.1 <- p.1 + theme(
               axis.ticks.x=element_blank(), 
               axis.line.x = element_line(colour = 'white', size=0.5, linetype='solid'),
               axis.line.y = element_line(colour = 'white', size=0.5, linetype='solid'),
               legend.position = "bottom", 
               legend.title=element_blank())
p.1 <- p.1 + labs(subtitle="10-FOLD 10-VALIDATION ON 91 IMMUNE PHENOTYPES",
              y="IMMUNE PHENOTYPES",
              x="DENSITY",
              title="DENSITY PREDICTION POWER",
              caption = "MATTHIJS KNIGGE")
p.1
ggsave(plot = p.1, filename = "Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/overview/predictiont.density.png", height = 15, width = 10, dpi = 300)



p.2 <- ggplot(data = cor.molt, aes(x = cytokines, fill = stimulation, color = stimulation))
p.2 <- p.2 + facet_grid(stimulation~., scales = "free", space = "free")
# p.2 <- p.2 + geom_jitter(col="grey", alpha=.7, show.legend = FALSE)
# p.2 <- p.2 + geom_density_ridges(alpha = 0.6, scale = .6)
# p.2 <- p.2 + geom_hline(yintercept = .05, size=.3, linetype="solid", color="red")
p.2 <- p.2 + geom_bar(stat="count", width = .4)
p.2 <- p.2 + theme(
  axis.ticks.x=element_blank(), 
  axis.line.x = element_line(colour = 'white', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'white', size=0.5, linetype='solid'),
  legend.position = "bottom", 
  legend.title=element_blank())
p.2 <- p.2 + labs(subtitle="10-FOLD 10-VALIDATION ON 91 IMMUNE PHENOTYPES",
              x="IMMUNE PHENOTYPES",
              y="AMOUNT OF GENES",
              title="AMOUNT OF GENES USED FOR MODEL BUILDING",
              caption = "MATTHIJS KNIGGE")
p.2 <- p.2 + coord_flip()
ggsave(plot = p.2, filename = "Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/overview/prediction.correlated.genes.png", width = 10, height = 15, dpi = 300)



c.mean$stimulation <- factor(c.mean$stimulation, levels = c("IFNy", "IL17", "IL1b", "IL22", "IL6", "TNFA"))

head(coef.molt)

p.3 <- ggplot(data = c.mean, aes(x = cytokine, y=gene.amount, fill = stimulation, color = stimulation))
p.3 <- p.3 + facet_grid(stimulation~., scales = "free", space = "free")
# p.3 <- p.3 + geom_jitter(col="grey", alpha=.7, show.legend = FALSE)
# p.3 <- p.3 + geom_density_ridges(alpha = 0.6, scale = .6)
# p.3 <- p.3 + geom_hline(yintercept = .05, size=.3, linetype="solid", color="red")
p.3 <- p.3 + geom_bar(stat="identity", width = .4)
p.3 <- p.3 + theme(
  axis.ticks.x=element_blank(), 
  axis.line.x = element_line(colour = 'white', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'white', size=0.5, linetype='solid'),
  legend.position = "bottom", 
  legend.title=element_blank())
p.3 <- p.3 + labs(subtitle="10-FOLD 10-VALIDATION ON 91 IMMUNE PHENOTYPES",
                  x="IMMUNE PHENOTYPES",
                  y="AMOUNT OF GENES",
                  title="AVERAGE AMOUNT GENES IN PREDICTED MODELS",
                  caption = "MATTHIJS KNIGGE")
p.3 <- p.3 + coord_flip()
ggsave(plot = p.3, filename = "Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/overview/predicted.genes.used.png", width = 10, height = 15, dpi = 300)



library(GGally)
pm <- ggmatrix(
  list(p, p.1, p.2, p.3),
  nrow = 1, ncol = 4,
  xAxisLabels = c("SPEARMEN COEFFICIENT", "DENSITY", "AMOUNT OF SIGNFICANT GENES", "AVERAGE AMOUNT OF GENES USED"),
  title = "10-FOLD 10-VALIDATION ON 91 IMMUNE PHENOTYPES"
)
ggsave(plot = pm, filename = "Bioinformatics/test.results/10.imputation.10.fold.immune.phenotypes/overview/overview.png", height = 15, width = 20, dpi = 300)













