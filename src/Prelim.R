### PRELIMINARY ANALYSIS

# Load Libraries
library(qtl)
library(ggplot2)
library(qtlnet)

# INSTALL QTLNET
# source("https://bioconductor.org/biocLite.R")
# biocLite("graph")
# biocLite("pcalg")

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
#phenos["Fat.wt"] <- "Fat.wt";
#phenos["LDL"] <- "LDL";
#phenos["HDL"] <- "HDL";
#phenos["CHOL"] <- "CHOL";
#phenos["adipose.turnover"] <- "adipose.turnover";
#phenos["liver.TG"] <- "liver.TG";
#phenos["Liver.wt"] <- "Liver.wt";
#phenos["INS.8wk"] <- "INS.8wk";
#phenos["TNF.alpha"] <- "TNF.alpha"; 
#phenos["GLU.4wk"] <- "GLU.4wk";
#phenos["GLU.6wk"] <- "GLU.6wk";
#phenos["GLU.8wk"] <- "GLU.8wk";
f2g$adipose <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], cbind(adipose.rz[,geneExps], phenotypes.rz[,phenos]))
f2g$liver <- cbind(f2g$pheno[,c("MouseNum", "Sex", "pgm")], cbind(liver.rz[,geneExps], phenotypes.rz[,phenos]))
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum", "Sex", "pgm")], cbind(adipose.rz[,geneExps], phenotypes.rz[,phenos]))
names(f2g$adipose) <- c(c("MouseNum", "Sex", "pgm"), names(geneExps), names(phenos))
names(f2g$liver) <- c(c("MouseNum", "Sex", "pgm"), names(geneExps), names(phenos))
names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), names(geneExps), names(phenos))
names(f2g$adipose)
names(f2g$liver)
names(f2g$pheno)

graphics.off()

# Pair-Wise plot of all the information
x11()
pairs(cbind(adipose.rz[,geneExps], (phenotypes[,phenos])), main="Adipose", labels=c(names(geneExps), names(phenos)), upper.panel=panel.cor, diag.panel=panel.hist)
x11()
pairs(cbind(liver.rz[,geneExps], (phenotypes[,phenos])), main="Liver", labels=c(names(geneExps), names(phenos)), upper.panel=panel.cor, diag.panel=panel.hist)

# Genome Scans

### LAZY SCANS! RUN ONCE AND SAVE!!! ###
f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

sex <- as.numeric(f2g$pheno$Sex) - 1

# CHANGE THIS NEED ALL PHENO COLUMNS!!
#f2g.fbn1.perm <- scanone(f2g, pheno.col=4:14, addcovar=sex, method="hk", n.perm=1000, perm.Xsp=TRUE)
#f2g.fbn1.perm2 <- scanone(f2g, pheno.col=4:6, addcovar=sex, method="hk", n.perm=1000, perm.Xsp=TRUE)
#f2g.fbn1.perm3 <- scanone(f2g, pheno.col=4:length(f2g$pheno), addcovar=sex, method="hk", n.perm=1000, perm.Xsp=TRUE)
setwd("data")
#save(list="f2g.fbn1.perm", file="Prelim-f2g.fbn1.perm.RData")
#save(list="f2g.fbn1.perm2", file="Prelim-f2g.fbn1.perm2.RData")
#save(list="f2g.fbn1.perm3", file="Prelim-f2g.fbn1.perm3.RData")
load(file="Prelim-f2g.fbn1.perm.RData")
load(file="Prelim-f2g.fbn1.perm2.RData")
load(file="Prelim-f2g.fbn1.perm3.RData")
setwd("..")

# CHOOSE THE RIGHT PERMS!!
PERMS <- f2g.fbn1.perm3

f2g.scan1 <- scanone(f2g, pheno.col=4:length(f2g$pheno), addcovar=sex, method="hk")

graphics.off()
for (i in 4:length(f2g$pheno)) {
  col <- i - 3
  plot(f2g.scan1, lodcolumn=col, main=paste(names(f2g$pheno)[i], "Scanone Plot"))
  add.threshold(f2g.scan1, perms=PERMS[,col], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,col], alpha=0.63, lty="dashed", lwd=2, col="green")
}

graphics.off()
for (i in 5:length(f2g$pheno)) {
  col <- i - 3
  f2g.scanAdd <- scanone(f2g, pheno.col=4, addcovar=f2g$pheno[,i], method="hk")
  f2g.scanInt <- scanone(f2g, pheno.col=4, intcovar=f2g$pheno[,i], method="hk")
  #x11()
  #par(mfrow=c(3,1))
  #plot(f2g.scan1, lodcolumn=1, main="FBN1")
  #plot(f2g.scanAdd, lodcolumn=1, main=paste("FBN1 with Additive Covariate", names(f2g$pheno)[i]))
  #plot(f2g.scanInt, lodcolumn=1, main=paste("FBN1 with Interactive Covariate", names(f2g$pheno)[i]))
  #par(mfrow=c(1,1))
  x11()
  plot(f2g.scan1, f2g.scanAdd, f2g.scanInt, lodcolumn=1, main=paste("FBN1 with Additive and Interactive Covariate", names(f2g$pheno)[i]))
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.63, lty="dashed", lwd=2, col="green")
}

### RUN QTL NET
FBN1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
SIRT1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Sirt1")]]
BMAL1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Arntl")]]
PIK3CG.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Pik3cg")]]
NR1H3.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nr1h3")]]
NCOR2.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ncor2")]]
SREBF1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Srebf1")]]
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("CHOL", "GLU.8wk", "TRIG.8wk", "TG.homogenate", "Liver.wt", "Leptin")],FBN1.adipose, SIRT1.adipose, BMAL1.adipose, PIK3CG.adipose, NR1H3.adipose, NCOR2.adipose, SREBF1.adipose)
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[c("CHOL", "GLU.8wk", "TRIG.8wk", "TG.homogenate", "Liver.wt", "Leptin")],FBN1.adipose, BMAL1.adipose)
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],FBN1.adipose, SIRT1.adipose, BMAL1.adipose, PIK3CG.adipose, NR1H3.adipose, NCOR2.adipose, SREBF1.adipose)
names(f2g$pheno)
