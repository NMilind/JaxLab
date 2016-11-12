### PRELIMINARY ANALYSIS

# Load Libraries
library(qtl)
library(ggplot2)

# Clear environmental variables
rm(list=ls())

# Load important functions
setwd("src")
source("important_func.R")
setwd("..")

# Load the BTBR Dataset
setwd("data")
load("BTBR.clean.data.Rdata")
setwd("..")
ls()

# The summary of the cross
summary(f2g)
# The names of the phenotype data currently being tested
names(f2g$pheno)
# Load some expression data into the phenotype list
# Analysis of SIRT1, SIRT6, CLOCK, and BMAL1
geneExps <- c()
phenos <- c()
geneExps["FBN1"] <- "1011565";
geneExps["SIRT1"] <- "519627";
geneExps["BMAL1"] <- "517966";
geneExps["PIK3CG"] <- "10002918744";
geneExps["NR1H3"] <- "505407";
geneExps["NCOR2"] <- "503860";
geneExps["SREBF1"] <- "10002908270";
geneExps["PER1"] <- "10002906335";
geneExps["PER2"] <- "501035";
geneExps["PER3"] <- "501565";
phenos["Fat.wt"] <- "Fat.wt";
phenos["GLU.4wk"] <- "GLU.4wk";
phenos["GLU.6wk"] <- "GLU.6wk";
phenos["GLU.8wk"] <- "GLU.8wk";
f2g$adipose <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], cbind(adipose.rz[,geneExps], phenotypes.rz[,phenos]))
f2g$liver <- cbind(f2g$pheno[,c("MouseNum", "Sex", "pgm")], cbind(liver.rz[,geneExps], phenotypes.rz[,phenos]))
names(f2g$adipose) <- c(c("MouseNum", "Sex", "pgm"), names(geneExps), names(phenos))
names(f2g$liver) <- c(c("MouseNum", "Sex", "pgm"), names(geneExps), names(phenos))
names(f2g$adipose)
names(f2g$liver)

graphics.off()

# Pair-Wise plot of all the information
x11()
pairs(cbind(adipose.rz[,geneExps], log(phenotypes[,phenos])), labels=c(names(geneExps), names(phenos)), upper.panel=panel.cor, diag.panel=panel.hist)
x11()
pairs(cbind(liver.rz[,geneExps], log(phenotypes[,phenos])), labels=c(names(geneExps), names(phenos)), upper.panel=panel.cor, diag.panel=panel.hist)
