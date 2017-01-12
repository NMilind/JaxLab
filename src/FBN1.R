#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

#################################################
## ENVIRONMENTAL SETUP                         ##
#################################################

source("src/configuration.R")

#################################################
## CROSS SETUP (ADIPOSE)                       ##
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
## SCANS                                       ##
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
graphics.open()
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
## PLOTS                                       ##
#################################################

graphics.off()
# Compare phenotypes to genotypes at identified QTL
graphics.open()
par(mfrow=c(2,2))
plot.pxg(f2g, find.marker(f2g, chr=1, pos=76.64), pheno.col=4, main="FBN1 QTL Chr1@76.64 cM")
plot.pxg(f2g, find.marker(f2g, chr=3, pos=49.12), pheno.col=4, main="FBN1 QTL Chr3@49.12 cM")
plot.pxg(f2g, find.marker(f2g, chr=7, pos=6.14), pheno.col=4, main="FBN1 QTL Chr7@6.14 cM")
plot.pxg(f2g, find.marker(f2g, chr="X", pos=8.60), pheno.col=4, main="FBN1 QTL ChrX@8.60 cM")
par(mfrow=c(1,1))

# Effect plots of the QTL
graphics.open()
par(mfrow=c(2,2))
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=1, pos=76.64), main="FBN1 QTL Chr1@76.64 cM")
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=3, pos=49.12), main="FBN1 QTL Chr3@49.12 cM")
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=7, pos=6.14), main="FBN1 QTL Chr7@6.14 cM")
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr="X", pos=8.60), main="FBN1 QTL ChrX@8.60 cM")
par(mfrow=c(1,1))

# Effect plots by sex
graphics.open()
effectplot(f2g, pheno.col=4, mname1="Sex", mark1=f2g$pheno$Sex, mname2=find.marker(f2g, chr=3, pos=49.12), main="FBN1 QTL Chr3@49.12 cM")
effectplot(f2g, pheno.col=4, mname1="Sex", mark1=f2g$pheno$Sex, mname2=find.marker(f2g, chr=7, pos=6.14), main="FBN1 QTL Chr7@6.14 cM")

# Obtain confidence intervals for QTL LOD peaks
graphics.open()
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
## MULTIPLE QTL ANALYSIS                       ##
#################################################

graphics.off()

# Scan for FBN1 in the cross with QTL Chr3 as additive covariate
FBN1.adipose.scan2 <- scanone(f2g, pheno.col=4, addcovar=f2g$geno$'3'$data[,find.marker(f2g, chr=3, pos=49.12)], method="hk")

# Plot the LOD scan with thresholds from the permutations
graphics.open()
plot(FBN1.adipose.scan2, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr7
# Interesting Peaks: ChrX

# Scan for FBN1 in the cross with QTL Chr7 as additive covariate
FBN1.adipose.scan3 <- scanone(f2g, pheno.col=4, addcovar=f2g$geno$'7'$data[,find.marker(f2g, chr=7, pos=6.14)], method="hk")

# Plot the LOD scan with thresholds from the permutations
graphics.open()
plot(FBN1.adipose.scan3, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr3

# QTL Effect Plots
graphics.off()
graphics.open()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=1, pos=76.64), mname1=find.marker(f2g, chr=3, pos=49.12), main="Chr1@76.64 cM x Chr3@49.12 cM")
graphics.open()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=1, pos=76.64), mname1=find.marker(f2g, chr=7, pos=6.14), main="Chr1@76.64 cM x Chr7@6.14 cM")
graphics.open()
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=1, pos=76.64), mname2=find.marker(f2g, chr="X", pos=8.60), main="Chr1@76.64 cM X ChrX@8.60 cM")
graphics.open()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=3, pos=49.12), mname1=find.marker(f2g, chr=7, pos=6.14), main="Chr3@49.12 cM x Chr7@6.14 cM")
graphics.open()
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=3, pos=49.12), mname2=find.marker(f2g, chr="X", pos=8.60), main="Chr3@49.12 cM x ChrX@8.60 cM")
graphics.open()
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
graphics.open()
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

#################################################
## ENVIRONMENTAL SETUP                         ##
#################################################

#source("src/configuration.R")

#################################################
## CROSS SETUP (HEPATIC)                       ##
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
FBN1.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
f2g$pheno <- cbind(f2g$pheno[,names(f2g$pheno)], FBN1.liver)
#  Check to make sure FBN1 is added
names(f2g$pheno)

#################################################
## SCANS                                       ##
#################################################

f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

sex <- as.numeric(f2g$pheno$Sex) - 1

# Permute through the cross to create baseline LOD score levels
#FBN1.liver.perms <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk", n.perm=1000, perm.Xsp=TRUE)
#save(file="data/FBN1-FBN1.liver.perms.RData", FBN1.liver.perms)
load(file="data/FBN1-FBN1.liver.perms.RData")

# Scan for FBN1 in the cross
FBN1.liver.scan1 <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk")

# Plot the LOD scan with thresholds from the permutations
graphics.open()
plot(FBN1.liver.scan1, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Interesting Peaks: Chr8, Chr12, and Chr15

# Tabulate and print LOD peaks where alpha=0.05
summary(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.05, format="tabByCol", ci.function="lodint")

# Tabulate and print LOD peaks where alpha=0.63
summary(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.63, format="tabByCol", ci.function="lodint")
# OUTPUT
#            chr   pos ci.low ci.high  lod
# rs13479718   8 28.03  17.92    54.9 3.82
# rs13481295  12  4.86   4.77    34.6 2.52
# rs13482700  15 43.63  31.94    55.7 3.24

#################################################
## MULTIPLE QTL ANALYSIS                       ##
#################################################

# Scan with Chr15 QTL as a covar to increase power on Chr8
FBN1.liver.scan2 <- scanone(f2g, pheno.col=4, addcovar=f2g$geno$'15'$data[,find.marker(f2g, chr=15, pos=43.63)], method="hk")

# Plot scan results
graphics.open()
plot(FBN1.liver.scan2, lodcolumn=1, main="Liver Tissue with Chr15 Covariate")
add.threshold(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr8
# Interesting Peaks: Chr12

# Tabulate and print LOD peaks where alpha=0.05
summary(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.05, format="tabByCol", ci.function="lodint")
# OUTPUT
#            chr pos ci.low ci.high  lod
# rs13479719   8  28   19.9    49.3 4.31

# Tabulate and print LOD peaks where alpha=0.63
summary(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.63, format="tabByCol", ci.function="lodint")
# OUTPUT
#            chr   pos ci.low ci.high  lod
# rs13479719   8 28.03  19.92    49.3 4.31
# rs13481295  12  4.86   4.77    33.5 3.10

# Compare phenotypes to genotypes at identified QTL
graphics.open()
plot.pxg(f2g, find.marker(f2g, chr=8, pos=28), pheno.col=4, main="FBN1 QTL Chr8@28.00 cM")

# Effect plots of the QTL
graphics.open()
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=8, pos=28), main="FBN1 QTL Chr8@28.00 cM")

# Effect plots by sex
graphics.open()
effectplot(f2g, pheno.col=4, mname1="Sex", mark1=f2g$pheno$Sex, mname2=find.marker(f2g, chr=8, pos=28), main="FBN1 QTL Chr8@28.00 cM")

# Obtain confidence intervals for QTL LOD peaks
graphics.open()
CI.Chr8 <- bayesint(FBN1.liver.scan2, chr=8, prob=0.95)
plot(FBN1.liver.scan2, chr=8, lodcolumn=1, main="Confidence Interval for Chr8")
lines(x=CI.Chr8[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
print(CI.Chr8[c(1,3),2])
# OUTPUT
# 20.92 cm to 48.92 cM

# Run a two-qtl scan over f2g
#gc()
#FBN1.liver.scan3 <- scantwo(f2g, pheno.col=4, addcovar=sex, method="hk")
#save(file="data/FBN1-FBN1.liver.scan3.RData", FBN1.liver.scan3)
load(file="data/FBN1-FBN1.liver.scan3.RData")

# Summary of results from scan
summary(FBN1.liver.scan3, c(9.1, 7.1, 6.3, 6.3, 3.3), "best")
#        pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
# c8:c15  27.9  43.8     7.83    4.01   0.561      27.9  43.8    7.27    3.45

graphics.off()

# Plot the results of the two-qtl scan
graphics.open()
plot(FBN1.liver.scan3)

# Make a QTL Model
FBN1.liver.model <- makeqtl(f2g, c("8","15"), c(27.9, 43.8))
FBN1.liver.modelFit <- fitqtl(f2g, pheno.col=4, qtl=FBN1.liver.model, formula=y~Q1+Q2)
summary(FBN1.liver.modelFit)
# ANOVA OUTPUT
# Variation accounted for: 6.232538%
# Model: y~Chr8+Chr15

# Two Additive QTL
# Chr8@27.9 cM 
# Chr15@43.8 cM

#################################################
## ENVIRONMENTAL SETUP                         ##
#################################################

#source("src/configuration.R")

#################################################
## CROSS SETUP (MUSCLE)                        ##
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
FBN1.gastroc <- gastroc.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
f2g$pheno <- cbind(f2g$pheno[,names(f2g$pheno)], FBN1.gastroc)
#  Check to make sure FBN1 is added
names(f2g$pheno)

#################################################
## SCANS                                       ##
#################################################

f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

sex <- as.numeric(f2g$pheno$Sex) - 1

# Permute through the cross to create baseline LOD score levels
#FBN1.gastroc.perms <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk", n.perm=1000, perm.Xsp=TRUE)
#save(file="data/FBN1-FBN1.gastroc.perms.RData", FBN1.gastroc.perms)
load(file="data/FBN1-FBN1.gastroc.perms.RData")

# Scan for FBN1 in the cross
FBN1.gastroc.scan1 <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk")

# Plot the LOD scan with thresholds from the permutations
graphics.open()
plot(FBN1.gastroc.scan1, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr9, Chr17
# Interesting Peaks: Chr2, Chr10, Chr16

# Tabulate and print LOD peaks where alpha=0.05
summary(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.05, format="tabByCol", ci.function="lodint")
# OUTPUT
#            chr  pos ci.low ci.high  lod
# rs13480263   9 40.1  20.95    52.4 4.39
# rs6358703   17 14.6   8.15    19.0 5.00

# Tabulate and print LOD peaks where alpha=0.63
summary(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.63, format="tabByCol", ci.function="lodint")
# OUTPUT
#            chr   pos ci.low ci.high  lod
# rs13476470   2 25.52   5.00    33.2 3.82
# rs13480263   9 40.11  20.95    52.4 4.39
# rs13480685  10 50.56  17.06    68.2 2.68
# rs4152386   16  2.35   2.35    20.2 2.71
# rs6358703   17 14.61   8.15    19.0 5.00

#################################################
## PLOTS                                       ##
#################################################

graphics.off()
# Compare phenotypes to genotypes at identified QTL
graphics.open()
par(mfrow=c(1,2))
plot.pxg(f2g, find.marker(f2g, chr=9, pos=40.10), pheno.col=4, main="FBN1 QTL Chr9@40.10 cM")
plot.pxg(f2g, find.marker(f2g, chr=17, pos=14.60), pheno.col=4, main="FBN1 QTL Chr17@14.60 cM")
par(mfrow=c(1,1))

# Effect plots of the QTL
graphics.open()
par(mfrow=c(1,2))
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=9, pos=40.10), main="FBN1 QTL Chr9@40.10 cM")
effectplot(f2g, pheno.col=4, mname1=find.marker(f2g, chr=17, pos=14.60), main="FBN1 QTL Chr17@14.60 cM")
par(mfrow=c(1,1))

# Effect plots by sex
graphics.open()
effectplot(f2g, pheno.col=4, mname1="Sex", mark1=f2g$pheno$Sex, mname2=find.marker(f2g, chr=9, pos=40.10), main="FBN1 QTL Chr9@40.10 cM")
effectplot(f2g, pheno.col=4, mname1="Sex", mark1=f2g$pheno$Sex, mname2=find.marker(f2g, chr=17, pos=14.60), main="FBN1 QTL Chr17@14.60 cM")

# Obtain confidence intervals for QTL LOD peaks
graphics.open()
par(mfrow=c(1,2))
CI.Chr9 <- bayesint(FBN1.gastroc.scan1, chr=9, prob=0.95)
plot(FBN1.gastroc.scan1, chr=9, lodcolumn=1, main="Confidence Interval for Chr9")
lines(x=CI.Chr9[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
print(CI.Chr9[c(1,3),2])
# From 21.45800 cM to 51.56976 cM
CI.Chr17 <- bayesint(FBN1.gastroc.scan1, chr=17, prob=0.95)
plot(FBN1.gastroc.scan1, chr=17, lodcolumn=1, main="Confidence Interval for Chr17")
lines(x=CI.Chr17[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
print(CI.Chr17[c(1,3),2])
# From 8.518018 cM to 17.332000 cM
par(mfrow=c(1,1))

#################################################
## MULTIPLE QTL ANALYSIS                       ##
#################################################

graphics.off()

# Scan for FBN1 in the cross with QTL Chr17 as additive covariate
FBN1.gastroc.scan2 <- scanone(f2g, pheno.col=4, addcovar=f2g$geno$'17'$data[,find.marker(f2g, chr=17, pos=14.60)], method="hk")

# Plot the LOD scan with thresholds from the permutations
graphics.open()
plot(FBN1.gastroc.scan2, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.gastroc.scan2, perms=FBN1.gastroc.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.gastroc.scan2, perms=FBN1.gastroc.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.gastroc.scan2, perms=FBN1.gastroc.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr9

# Scan for FBN1 in the cross with QTL Chr9 as additive covariate
FBN1.gastroc.scan3 <- scanone(f2g, pheno.col=4, addcovar=f2g$geno$'9'$data[,find.marker(f2g, chr=9, pos=40.10)], method="hk")

# Plot the LOD scan with thresholds from the permutations
graphics.open()
plot(FBN1.gastroc.scan3, lodcolumn=1, main="FBN1 LOD for Variation across Mouse Genome")
add.threshold(FBN1.gastroc.scan3, perms=FBN1.gastroc.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.gastroc.scan3, perms=FBN1.gastroc.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.gastroc.scan3, perms=FBN1.gastroc.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# Major Peaks: Chr17, Chr2

# QTL Effect Plots
graphics.off()
graphics.open()
effectplot(f2g, pheno.col=4, mname2=find.marker(f2g, chr=9, pos=40.10), mname1=find.marker(f2g, chr=17, pos=14.60), main="Chr9@40.10 cM x Chr17@14.60 cM")

# Run a two-qtl scan over f2g
#gc()
#FBN1.gastroc.scan4 <- scantwo(f2g, pheno.col=4, addcovar=sex, method="hk")
#save(file="data/FBN1-FBN1.gastroc.scan4.RData", FBN1.gastroc.scan4)
#load(file="data/FBN1-FBN1.gastroc.scan4.RData")

# Summary of results from scan
#summary(FBN1.adipose.scan4, c(9.1, 7.1, 6.3, 6.3, 3.3), "best")
#       pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
# c3:c7  62.5  8.14     9.66    5.35    1.32      49.5  6.14    8.34    4.03

#graphics.off()

# Plot the results of the two-qtl scan
#graphics.open()
#plot(FBN1.adipose.scan4)

# Make a QTL Model
#FBN1.adipose.model <- makeqtl(f2g, c("3","7"), c(49.5, 6.14))
#FBN1.adipose.modelFit <- fitqtl(f2g, pheno.col=4, qtl=FBN1.adipose.model, formula=y~Q1+Q2)
#summary(FBN1.adipose.modelFit)
# ANOVA OUTPUT
# Variation accounted for: 7.384149%
# Model: y~Chr3+Chr7

# Two Additive QTL
# Chr3@49.5 cM 
# Chr7@6.14 cM

#################################################
## TANDEM GRAPHS                               ##
#################################################

graphics.off()

par(mfrow=c(3,1))

plot(FBN1.gastroc.scan1, lodcolumn=1, main="Muscle Tissue")
add.threshold(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.gastroc.scan1, perms=FBN1.gastroc.perms, alpha=0.63, lty="dashed", lwd=1, col="red")

plot(FBN1.adipose.scan1, lodcolumn=1, main="Adipose Tissue")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")

plot(FBN1.adipose.scan2, lodcolumn=1, main="Adipose Tissue with Chr3 Covariate")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan2, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")

plot(FBN1.adipose.scan3, lodcolumn=1, main="Adipose Tissue with Chr7 Covariate")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan3, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")

plot(FBN1.liver.scan1, lodcolumn=1, main="Liver Tissue")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.63, lty="dashed", lwd=1, col="red")

plot(FBN1.liver.scan2, lodcolumn=1, main="Liver Tissue with Chr15 Covariate")
add.threshold(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.liver.scan2, perms=FBN1.liver.perms, alpha=0.63, lty="dashed", lwd=1, col="red")

par(mfrow=c(1,1))

par(mfrow=c(2,1))
plot(FBN1.adipose.scan1, lodcolumn=1, main="Adipose Tissue")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.adipose.scan1, perms=FBN1.adipose.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
plot(FBN1.liver.scan1, lodcolumn=1, main="Liver Tissue")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.05, lty="dashed", lwd=1, col="green")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
add.threshold(FBN1.liver.scan1, perms=FBN1.liver.perms, alpha=0.63, lty="dashed", lwd=1, col="red")
par(mfrow=c(1,1))
