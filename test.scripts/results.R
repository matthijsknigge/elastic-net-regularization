# library
library(elastic.net.regularization)
library(data.table)
library(GGally)


# CYTOKINES
cytokines <- fread("~/Bioinformatics/results/10.imputation.10.fold/summary/cytokines.coef.txt")
# CELL COUNTS
cell.counts <- fread("~/Bioinformatics/results/10.imputation.10.fold/summary/cell.counts.coef.txt")
# PLATELET COUNTS
platelet.counts <- fread("~/Bioinformatics/results/10.imputation.10.fold/summary/platelet.counts.coef.txt")
# HORMONE LEVELS
hormone.levels <- fread("~/Bioinformatics/results/10.imputation.10.fold/summary/hormone.levels.coef.txt")
# PROTEIN LEVELS
protein.levels <- fread("~/Bioinformatics/results/10.imputation.10.fold/summary/protein.levels.coef.txt")
# IMMUNOGLOBULIN
immunoglobulin <- fread("~/Bioinformatics/results/10.imputation.10.fold/summary/immunoglobulin.levels.coef.txt")
# METABOLITES
metabolites <- fread("~/Bioinformatics/results/10.imputation.10.fold/summary/metabolites.coef.txt")
# summary
summary <- rbind(cytokines, cell.counts, platelet.counts, hormone.levels, protein.levels, immunoglobulin, metabolites)

# plotting
for(name in c("cytokines", "cell.counts", "platelet.counts", "hormone.levels", "protein.levels", "immunoglobulin", "metabolites")){
  data <- get(name)
  data <- glmnet.create.mean.labels(data = data)$d
  p <- glmnet.box.plot(coef = data$coef,
                       group = data$group,
                       identifier = data$identifier,
                       title = name)$p

  h <- length(unique(data$identifier))*.15
  if(length(unique(data$identifier)) < 20){
    h <- 3
  }
  print(h)
  ggsave(filename = paste0("~/Bioinformatics/results/10.imputation.10.fold/summary/", name, ".boxplot.png"), dpi = 300,
         height = h, width = 10, plot = p)

}

names <- c("cell.counts", "cytokines", "hormone.levels", "immunoglobulin.levels", "metabolites", "platelet.counts", "protein.levels")

# data
for(name in names){

  # box
  path <- "~/Bioinformatics/results/10.imputation.10.fold/summary/"
  data <- fread(paste0(path, name, "/", name, ".coef.melt.txt"))
  data <- glmnet.create.mean.labels(data = data)$d
  p.1 <- glmnet.box.plot(coef = data$coef, group = data$group, identifier = data$identifier, title = name)$p

  # density
  data <- fread(paste0(path, name, "/", name, ".coef.melt.txt"))
  p.2 <- glmnet.density.plot(coef = data$coef, group = data$group, identifier = data$identifier, title = name)$p

  # bar
  data <- fread(paste0(path, name, "/", name, ".coef.txt"))
  p.3 <- glmnet.bar.plot(identifier = data$identifier, group = data$group, gene.amount.before = data$gene.amount, gene.amount.after = data$gene.amount.significant, title = name)$p

  # grid
  pm <- ggmatrix(
    list(p.1, p.2, p.3),
    nrow = 1, ncol = 3,
    xAxisLabels = c("SPEARMEN COEFFICIENT", "DENSITY", "AMOUNT OF GENES"),
    title = name
  )

  h <- length(unique(data$identifier))*.15

  if(length(unique(data$identifier)) < 20){
    h <- 3
  }

  ggsave(plot = pm, filename = paste0("~/Bioinformatics/results/10.imputation.10.fold/summary/", name, ".png"),
         dpi = 300, height = h, width = 10)

}

