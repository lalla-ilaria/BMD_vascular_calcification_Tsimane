#libraries
library(rethinking)
library(tidyverse)

#plotting parameters
malecol <- "turquoise4"
femalecol <- "darkorange"
allcol <- "purple4"
directcol <- "skyblue2"
partialcol1 <- "cornflowerblue"
partialcol2 <- "dodgerblue3"
totalcol <- "blue4"
cols <- c(femalecol, malecol)
seqalong <- seq(-4,4, length.out = 100)
xlim <- c(-2, 2)
ylim <- c(-1, 1)
mar <- c(4, 3, 0.5, 0.5)
counterfactual_age <- 60
gq <- list() #generated quantities

source("0_functions.R")
source("1_data_process.R")
source("2_analysis.R")
source("3_outcomes.R")