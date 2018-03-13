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













