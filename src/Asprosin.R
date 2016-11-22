#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

#################################################
## 1: ENVIRONMENTAL SETUP                      ##
#################################################

# Clear environmental variables
rm(list=ls())

# Import required libraries
library(qtl)
library(ggplot2)

# Set working directory
setwd("~/Desktop/JaxLab")

# Import generic functions
source("src/important_func.R")

# Import BTBR data
load(file="data/BTBR.clean.data.Rdata")
