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

#################################################
## 2: CROSS SETUP (ADIPOSE)                    ##
#################################################

# Print a summary of the BTBR cross
summary(f2g)
# Print the names of the current phenotypes
names(f2g$pheno)
# Use only MouseNum, Sex, and pgm for analytics
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
# Check to make sure only three columns are conserved
names(f2g$pheno)
# Add FBN1 as a phenotype
FBN1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
f2g$pheno <- cbind(f2g$pheno[,names(f2g$pheno)], FBN1.adipose)
#  Check to make sure FBN1 is added
names(f2g$pheno)

#################################################
## 3: SCANS                                    ##
#################################################

f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

sex <- as.numeric(f2g$pheno$Sex) - 1

# Permute through the cross to create baseline LOD score levels
#FBN1.adipose.perms <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk", n.perm=1000, perm.Xsp=TRUE)
#save(file="data/FBN1-FBN1.adipose.perms.RData", FBN1.adipose.perms)
load(file="data/FBN1-FBN1.adipose.perms.RData")

# Scan for FBN1 in the cross
FBN1.adipose.scan1 <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk")

# Plot the LOD scan with thresholds from the permutations
x11()
plot(FBN1.adipose.scan1, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr3, Chr7
# Interesting Peaks: Chr1, ChrX

# Tabulate and print LOD peaks where alpha=0.05
summary(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.05, format="tabByCol", ci.function="lodint")
# OUTPUT
#            chr   pos ci.low ci.high  lod
# rs13477374   3 49.12  44.55    65.1 4.40
# c7.loc4      7  6.14   2.14    15.9 4.07

# Tabulate and print LOD peaks where alpha=0.63
summary(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.63, format="tabByCol", ci.function="lodint")
# OUTPUT
#            chr   pos ci.low ci.high  lod
# c1.loc75     1 76.64   1.64    86.6 2.83
# rs13477374   3 49.12  44.55    65.1 4.40
# c7.loc4      7  6.14   2.14    15.9 4.07
# cX.loc5      X  8.60   3.60    28.3 2.17

#################################################
## 4: PLOTS                                    ##
#################################################

graphics.off()
# Compare phenotypes to genotypes at identified QTL
x11()
par(mfrow=c(2,2))
plot.pxg(f2g, find.marker(f2g, chr=1, pos=76.64), pheno.col=4, main="FBN1 QTL Chr1@76.64 cM")
plot.pxg(f2g, find.marker(f2g, chr=3, pos=49.12), pheno.col=4, main="FBN1 QTL Chr3@49.12 cM")
plot.pxg(f2g, find.marker(f2g, chr=7, pos=6.14), pheno.col=4, main="FBN1 QTL Chr7@6.14 cM")
plot.pxg(f2g, find.marker(f2g, chr="X", pos=8.60), pheno.col=4, main="FBN1 QTL ChrX@8.60 cM")
par(mfrow=c(1,1))

# Effect plots of the QTL
x11()
par(mfrow=c(2,2))
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=1, pos=76.64), main="FBN1 QTL Chr1@76.64 cM")
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=3, pos=49.12), main="FBN1 QTL Chr3@49.12 cM")
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=7, pos=6.14), main="FBN1 QTL Chr7@6.14 cM")
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr="X", pos=8.60), main="FBN1 QTL ChrX@8.60 cM")
par(mfrow=c(1,1))

# Obtain confidence intervals for QTL LOD peaks
x11()
par(mfrow=c(2,2))
CI.Chr1 <- bayesint(FBN1.adipose.scan1, chr=1, prob=0.95)
plot(FBN1.adipose.scan1, chr=1, lodcolumn=1, main="Confidence Interval for Chr1")
lines(x=CI.Chr1[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
print(CI.Chr1[c(1,3),2])
# From 1.64490 cM to 88.34743 cM
CI.Chr3 <- bayesint(FBN1.adipose.scan1, chr=3, prob=0.95)
plot(FBN1.adipose.scan1, chr=3, lodcolumn=1, main="Confidence Interval for Chr3")
lines(x=CI.Chr3[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
print(CI.Chr3[c(1,3),2])
# From 46.54620 cM to 63.65573 cM
CI.Chr7 <- bayesint(FBN1.adipose.scan1, chr=7, prob=0.95)
plot(FBN1.adipose.scan1, chr=7, lodcolumn=1, main="Confidence Interval for Chr7")
lines(x=CI.Chr7[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
print(CI.Chr7[c(1,3),2])
# From 2.1403 cM to 15.1403 cM
CI.ChrX <- bayesint(FBN1.adipose.scan1, chr="X", prob=0.95)
plot(FBN1.adipose.scan1, chr="X", lodcolumn=1, main="Confidence Interval for ChrX")
lines(x=CI.ChrX[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
print(CI.ChrX[c(1,3),2])
# From 3.60400 cM to 26.20655 cM
par(mfrow=c(1,1))

#################################################
## 5: MULTIPLE QTL ANALYSIS                    ##
#################################################

graphics.off()

# Scan for FBN1 in the cross with QTL Chr3 as additive covariate
FBN1.adipose.scan2 <- scanone(f2g, pheno.col=4, addcovar=f2g$geno$'3'$data[,find.marker(f2g, chr=3, pos=49.12)], method="hk")

# Plot the LOD scan with thresholds from the permutations
x11()
plot(FBN1.adipose.scan2, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr7
# Interesting Peaks: ChrX

# Scan for FBN1 in the cross with QTL Chr7 as additive covariate
FBN1.adipose.scan3 <- scanone(f2g, pheno.col=4, addcovar=f2g$geno$'7'$data[,find.marker(f2g, chr=7, pos=6.14)], method="hk")

# Plot the LOD scan with thresholds from the permutations
x11()
plot(FBN1.adipose.scan3, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr3

# QTL Effect Plots
graphics.off()
x11()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=1, pos=76.64), mname1=find.marker(f2g, chr=3, pos=49.12), main="Chr1@76.64 cM x Chr3@49.12 cM")
x11()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=1, pos=76.64), mname1=find.marker(f2g, chr=7, pos=6.14), main="Chr1@76.64 cM x Chr7@6.14 cM")
x11()
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=1, pos=76.64), mname2=find.marker(f2g, chr="X", pos=8.60), main="Chr1@76.64 cM X ChrX@8.60 cM")
x11()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=3, pos=49.12), mname1=find.marker(f2g, chr=7, pos=6.14), main="Chr3@49.12 cM x Chr7@6.14 cM")
x11()
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=3, pos=49.12), mname2=find.marker(f2g, chr="X", pos=8.60), main="Chr3@49.12 cM x ChrX@8.60 cM")
x11()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=7, pos=6.14), mname1=find.marker(f2g, chr="X", pos=8.60), main="Chr7@6.14 cM x ChrX@8.60 cM")

# Run a two-qtl scan over f2g
#gc()
#FBN1.adipose.scan4 <- scantwo(f2g, pheno.col=4, addcovar=sex, method="hk")
#save(file="data/FBN1-FBN1.adipose.scan4.RData", FBN1.adipose.scan4)
load(file="data/FBN1-FBN1.adipose.scan4.RData")

# Summary of results from scan
summary(FBN1.adipose.scan4, c(9.1, 7.1, 6.3, 6.3, 3.3), "best")
#       pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
# c3:c7  62.5  8.14     9.66    5.35    1.32      49.5  6.14    8.34    4.03

graphics.off()

# Plot the results of the two-qtl scan
x11()
plot(FBN1.adipose.scan4)

# Make a QTL Model
FBN1.adipose.model <- makeqtl(f2g, c("3","7"), c(49.5, 6.14))
FBN1.adipose.modelFit <- fitqtl(f2g, pheno.col=4, qtl=FBN1.adipose.model, formula=y~Q1+Q2)
summary(FBN1.adipose.modelFit)
# ANOVA OUTPUT
# Variation accounted for: 7.384149%
# Model: y~Chr3+Chr7

# Two Additive QTL
# Chr3@49.5 cM 
# Chr7@6.14 cM