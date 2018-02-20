# DIMENSIONS

# N = sample size
# K = number of genetic variants

# INDIVIDUAL LEVEL DATA
# g = genetic variant(s), matrix dimension N x K
# x = risk factor/exposure, vector length n
# y = outcome, vector length n

# SUMMARIZED DATA
# Bx = genetic association with exposure, vector length K
# Bx.se = standard errors of genetic associations with exposure
# By = genetic association with outcome, vector length K
# By.se = standard errors of genetic associations with outcome


getwd()
setwd("/home/mknigge/Bioinformatics/interplay.gut.microbiome.immune.system/")
list.dirs()

library(devtools)
library(roxygen2)
create("cross.validation")
