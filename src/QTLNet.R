#############################################
#Jackson Laboratory Summer Student Program
# Example qtlnet code
#Attie BTBR eQTL data
##############################################
#
#set working directory
rm(list=ls())
directory <- "~/Desktop/JaxLab"
setwd(directory)
# load libraries
library(qtl)
library(ggplot2)
# load myfunctions.R script
setwd("src")
source("important_func.R")
setwd("..")
# load clean BTBR RData file
setwd("data")
load("BTBR.clean.data.Rdata")
setwd("..")
ls()
#
#define variables
FBN1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
SIRT1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Sirt1")]]
BMAL1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Arntl")]]
PIK3CG.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Pik3cg")]]
NR1H3.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nr1h3")]]
NCOR2.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ncor2")]]
SREBF1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Srebf1")]]
TGFB1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Tgfb1")]]

FBN1.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
SIRT1.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Sirt1")]]
BMAL1.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Arntl")]]
PIK3CG.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Pik3cg")]]
NR1H3.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nr1h3")]]
NCOR2.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ncor2")]]
SREBF1.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Srebf1")]]
TGFB1.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Tgfb1")]]

f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("CHOL", "GLU.8wk", "TRIG.8wk", "TG.homogenate", "Liver.wt", "Leptin")],FBN1.adipose, SIRT1.adipose, BMAL1.adipose, PIK3CG.adipose, NR1H3.adipose, NCOR2.adipose, SREBF1.adipose)
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("CHOL", "GLU.8wk", "TRIG.8wk", "TG.homogenate", "Liver.wt", "Leptin")],FBN1.adipose, BMAL1.adipose)
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],FBN1.adipose, SIRT1.adipose, BMAL1.adipose, PIK3CG.adipose, NR1H3.adipose, NCOR2.adipose, SREBF1.adipose, TGFB1.adipose)
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],FBN1.liver, SIRT1.liver, BMAL1.liver, PIK3CG.liver, NR1H3.liver, NCOR2.liver, SREBF1.liver, TGFB1.liver)
names(f2g$pheno)

setwd("net6-data")

# Step 0:  Load qtlnet
library(qtlnet)
#qtlnet of Ppy, Pyy, Npy (Extended Fig. 1)

#define variables
#f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],Ppy,Pyy,Npy)

################################################################
# Go through the 5 Phases of QTLnet (Script provided by Ying Qi)
################################################################

###############################################################################
#STEP 1: Set up for network with nodes pheno.col and maximum parent max.parents
###############################################################################
# Phenotype identifer from cross object.
pheno.col <- 4:ncol(f2g$pheno)
#
# maximum number of parents is 3...
max.parents <- ncol(f2g$pheno)-4
#
size.qtlnet(pheno.col, max.parents)
# computes the number of possible scanones with nodes pheno.col and max parent 
# 4* (3C1 + 3C2 + 3C3) + 4*(1)
#
parents <- parents.qtlnet(pheno.col, max.parents) 
# Creates a list of all parent sets to be used as covariates of child phenotypes in scanone calculations
#
groups <- group.qtlnet(parents = parents, group.size = 8)
# groups parents into approx equal sized groups for parallel computations
#
save(f2g, pheno.col, max.parents, parents, groups, file = "Step1.RData", compress = TRUE)
pa <- summary(parents)
# Save files for later computational use
#
########################################################################
# STEP 2: run scanone on each phenotype conditioned on parent phenotypes
########################################################################
# precompute BIC scores using scanone on a "grid of computers"; parallel computing
N <- rep(NA,nrow(groups))
# Creates an empty matrix 
#
for (i in 1:nrow(groups))
  N[i] <- sum(pa[seq(groups[i,1], groups[i,2]),2])
# Calculates total number of networks to per group,
# Assign this sum to N
#
# for each row in groups, calculate BIC score
for(i in seq(nrow(groups))){
  ##TODO: intcov. Type ?bic.qtlnet for more info
  bic <- bic.qtlnet(f2g,
                    pheno.col,
                    threshold = 3.04,
                    max.parents = max.parents,
                    parents = parents[seq(groups[i,1], groups[i,2])],
                    verbose = FALSE)
  bic
  # gets the list containing all possible parents
  #
  #########################
  # STEP 3: Save BIC scores
  #########################
  save(bic, file = paste("bic", i, ".RData", sep = ""), compress = TRUE)
}
###########################
# STEP 4: MCMC, parallelize
###########################
load("Step1.RData")
# Read in saved BIC score and combine into one object
bic.group <- list()
bic.group
for(i in seq(nrow(groups))){
  load(paste("bic", i, ".RData", sep = ""))
  bic.group[[i]] <- bic
  cat("group =", i, "\n")
}
#
saved.scores <- bic.join(f2g, pheno.col = 4:ncol(f2g$pheno), bic.group, max.parents = ncol(f2g$pheno)-4)
saved.scores
# sample Markov chain (MCMC)
set.seed(54321)
#
n.runs <- 5
# sets number of parallel runs
#
# parallelize this
for (i in seq(n.runs))
{ 
  cat("run = ", i, "\n")
  ## Run MCMC with randomized initial network
  mcmc <- mcmc.qtlnet(f2g,
                      pheno.col=4:ncol(f2g$pheno),
                      threshold = 3.04,
                      thinning = 1,  #thinning rate
                      max.parents = max.parents,
                      saved.scores = saved.scores, #precomputed bic score
                      init.edges = NULL ) #edges chosen uniformly from 0 to number of possible edges  
  save(mcmc, file = paste("mcmc", i, ".RData", sep = ""), compress = TRUE)
}
#############################################
# STEP 5: Combine results for Post-processing
#############################################
outs.qtlnet <-list()
for (i in seq(n.runs)){
  load(paste("mcmc", i, ".RData", sep = ""))
  outs.qtlnet[[i]] <- mcmc
}
#
out1 <-c.qtlnet(outs.qtlnet) # catenates the outputs from 3 separate runs together
#
summary(out1)
plot(out1)

save(file="out1.RData", out1)

###################################
# Opening data from previous runs
###################################
setwd("..")
setwd("net5-data")
load(file="out1.RData")
summary(out1)
plot(out1)

