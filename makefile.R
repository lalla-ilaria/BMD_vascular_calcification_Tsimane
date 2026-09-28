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


warmup <- 500
iter <- 1000

source("0_functions.R")
source("1_bis_simulate_data_to_test_code.R")
#if original data are available
#source("1_data_process.R")
#to use simulated data
data_to_use <- dat_sim
#to use real data (if available)
#data_to_use <- all_dat
source("2_analysis.R")
source("3_outcomes.R")
