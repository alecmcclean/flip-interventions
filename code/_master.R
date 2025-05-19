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
       latex2exp, 
       SuperLearner,
       wooldridge)

options(stringsAsFactors = F)

source("1_wage_data_analysis.R")
source("3_simulations.R")