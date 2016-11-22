# JaxLab
--------
Create a `data` folder in the root folder of the repository and put the BTBR data in there for the scripts to work. Please do not forget to change the paths to the directories at the beginning of each script file.

# Configuration
---------------
Please create a file called `configuration.R` in the `src` folder. Add the following content to the file:
```
# Configuration file

# Clear environmental variables
rm(list=ls())

# Import required libraries
library(qtl)
library(ggplot2)

# Set working directory
setwd("YOUR WORKING DIRECTORY")

# Import generic functions
source("src/important_func.R")

# Import BTBR data
load(file="data/BTBR.clean.data.Rdata")
```