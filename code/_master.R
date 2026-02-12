###############################################################################
### Author: Alec McClean
### Purpose: Master script for longitudinal flip interventions
###############################################################################

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(magrittr,
       tidyverse,
       ggthemes,
       glmnet,
       ranger,
       rpart,
       latex2exp, 
       SuperLearner,
       wooldridge)

options(stringsAsFactors = F)

# Set working directory here
this_file <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(this_file))

source("1_wage_data_analysis.R")
source("3_simulations.R")