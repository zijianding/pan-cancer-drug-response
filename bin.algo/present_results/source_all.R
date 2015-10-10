###library###
library(pracma)
library(doParallel)
library(foreach)
no_cores = detectCores()
library(pracma)
library(ggplot2)

#functions
setwd("functions/")
source("auc_random.R")
source("ave_roc.R")
source("ensemble_roc.R")
source("error_calc.R")
source("find_name.R")
source("identify_class.R")
source("read_table.R")
source("multiplot.R")