#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

#################################################
## ENVIRONMENTAL SETUP                         ##
#################################################

source("src/configuration.R")

#################################################
## CROSS SETUP                                 ##
#################################################

# Print a summary of the BTBR cross
summary(f2g)
# Print the names of the current phenotypes
names(f2g$pheno)
# Use only MouseNum, Sex, and pgm for analytics
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
# Check to make sure only three columns are conserved
names(f2g$pheno)
# Add GNG8 as a phenotype
GNG8.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Gng8")]]
f2g$pheno <- cbind(f2g$pheno[,names(f2g$pheno)], GNG8.liver)
#  Check to make sure GNG8 is added
names(f2g$pheno)

#################################################
## SCANS                                       ##
#################################################

f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

sex <- as.numeric(f2g$pheno$Sex) - 1

# Permute through the cross to create baseline LOD score levels
#GNG8.liver.perms <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk", n.perm=1000, perm.Xsp=TRUE)
#save(file="data/GNG8-GNG8.liver.perms.RData", GNG8.liver.perms)
load(file="data/GNG8-GNG8.liver.perms.RData")

# Scan for FBN1 in the cross
GNG8.liver.scan1 <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk")

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