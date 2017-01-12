### PRELIMINARY ANALYSIS

source("src/configuration.R")

# INSTALL QTLNET
# source("https://bioconductor.org/biocLite.R")
# biocLite("graph")
# biocLite("pcalg")

# The summary of the cross
summary(f2g)
# The names of the phenotype data currently being tested
names(f2g$pheno)
# Load some expression data into the phenotype list
# Analysis of SIRT1, SIRT6, CLOCK, and BMAL1
geneExps <- c()
phenos <- c()
geneExps["FBN1"] <- "1011565"
geneExps["GPR61"] <- "10002932198"
geneExps["GPSM2"] <- "505911"
geneExps["FNDC7"] <- "10002917017"
geneExps["COL11A1"] <- "501818"
geneExps["S1PR1"] <- "10002919050"
geneExps["GPR88"] <- "519210"
geneExps["AGL"] <- "516642"
geneExps["CNN3"] <- "10002908968"
geneExps["MYOZ2"] <- "505596"
geneExps["PLA2G12A"] <- "592553"
geneExps["COL25A1"] <- "10004033701"
geneExps["CYP2U1"] <- "10002913136"
geneExps["EMCN"] <- "10002916274"
#geneExps["SIRT1"] <- "519627";
#geneExps["BMAL1"] <- "517966";
#geneExps["PIK3CG"] <- "10002918744";
#geneExps["NR1H3"] <- "505407";
#geneExps["NCOR2"] <- "503860";
#geneExps["SREBF1"] <- "10002908270";
#geneExps["TGFB1"] <- "507240";
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
f2g$adipose <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")], cbind(liver.rz[,geneExps], phenotypes.rz[,phenos]))
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
quartz()
pairs(cbind(adipose.rz[,geneExps], (phenotypes[,phenos])), main="Adipose", labels=c(names(geneExps), names(phenos)), upper.panel=panel.cor, diag.panel=panel.hist)
quartz()
pairs(cbind(liver.rz[,geneExps], (phenotypes[,phenos])), main="Liver", labels=c(names(geneExps), names(phenos)), upper.panel=panel.cor, diag.panel=panel.hist)

# Pair-Wise plots using proper log transforms
FBN1.adipose <- adipose.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
Pheno.Fat.wt <- phenotypes[,"Fat.wt"]
Pheno.LDL <- log(phenotypes[,"LDL"])
Pheno.CHOL <- log(phenotypes[,"CHOL"])
Pheno.INS.8wk <- log(phenotypes[,"INS.8wk"])
Pheno.GLU.4wk <- log(phenotypes[,"GLU.4wk"])
Pheno.Weight <- phenotypes[,"Weight"]
Pheno.Gastroc.wt <- phenotypes[,"Gastroc.wt"]
pairs(
  cbind(FBN1.adipose, Pheno.Fat.wt, Pheno.LDL, Pheno.INS.8wk, Pheno.GLU.4wk, Pheno.Weight, Pheno.Gastroc.wt), 
  main="FBN1", 
  labels=c("FBN1", "Fat.wt", "LDL", "INS.8wk", "GLU.4wk", "Weight", "Gastroc.wt"), 
  upper.panel=panel.cor, 
  diag.panel=panel.hist
)

FBN1.liver <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Sirt1")]]
Pheno.Liver.wt <- log(phenotypes[,"Liver.wt"])
Pheno.HDL <- log(phenotypes[,"HDL"])
Pheno.Leptin <- log(phenotypes[,"Leptin"])
Pheno.liver.TG <- log(phenotypes[,"liver.TG"])
pairs(
  cbind(FBN1.adipose, FBN1.liver, Pheno.Liver.wt, Pheno.HDL, Pheno.Leptin, Pheno.liver.TG), 
  main="FBN1", 
  labels=c("FBN1.adipose", "FBN1.liver", "Liver.wt", "HDL", "Leptin", "liver.TG"), 
  upper.panel=panel.cor, 
  diag.panel=panel.hist
)

FBN1.gastroc <- gastroc.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
Pheno.Gastroc.wt <- phenotypes[,"Gastroc.wt"]
Pheno.Fat.wt <- phenotypes[,"Fat.wt"]
Pheno.LDL <- log(phenotypes[,"LDL"])
Pheno.CHOL <- log(phenotypes[,"CHOL"])
Pheno.INS.8wk <- log(phenotypes[,"INS.8wk"])
Pheno.GLU.4wk <- log(phenotypes[,"GLU.4wk"])
Pheno.Weight <- phenotypes[,"Weight"]
pairs(
  cbind(FBN1.gastroc, Pheno.Gastroc.wt, Pheno.Fat.wt, Pheno.LDL, Pheno.CHOL, Pheno.INS.8wk, Pheno.GLU.4wk, Pheno.Weight),
  main="FBN1",
  labels=c("FBN1.gastroc", "Gastroc.wt", "Fat.wt", "LDL", "CHOL", "INS.8wk", "GLU.4wk", "Weight"),
  upper.panel=panel.cor,
  diag.panel=panel.hist
)

FBN1.adipose <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Fbn1")]]
SIRT1.adipose <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Sirt")]]
PIK3CG.adipose <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Pik3cg")]]
NR1H3.adipose <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Nr1h3")]]
NCOR2.adipose <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Ncor2")]]
SREBF1.adipose <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Srebf1")]]
TGFB.adipose <- liver.rz[,annot$a_gene_id[which(annot$gene_symbol=="Tgfb1")]]
pairs(
  cbind(FBN1.adipose, SIRT1.adipose, PIK3CG.adipose, NR1H3.adipose, NCOR2.adipose, SREBF1.adipose, TGFB.adipose),
  main="FBN1",
  labels=c("Fbn1", "Sirt1", "Pik3cg", "Nr1h3", "Ncor2", "Srebf1", "Tgfb1"),
  upper.panel=panel.cor,
  diag.panel=panel.hist
)

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
  #quartz()
  #par(mfrow=c(3,1))
  #plot(f2g.scan1, lodcolumn=1, main="FBN1")
  #plot(f2g.scanAdd, lodcolumn=1, main=paste("FBN1 with Additive Covariate", names(f2g$pheno)[i]))
  #plot(f2g.scanInt, lodcolumn=1, main=paste("FBN1 with Interactive Covariate", names(f2g$pheno)[i]))
  #par(mfrow=c(1,1))
  #quartz()
  par(mfrow=c(3,1))
  plot(f2g.scan1, lodcolumn=1, main="FBN1")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.63, lty="dashed", lwd=2, col="green")
  plot(f2g.scan1, f2g.scanAdd, lodcolumn=1, main=paste("FBN1 with Additive Covariate", names(f2g$pheno)[i]))
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.63, lty="dashed", lwd=2, col="green")
  plot(f2g.scan1, f2g.scanInt, lodcolumn=1, main=paste("FBN1 with Interactive Covariate", names(f2g$pheno)[i]))
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.63, lty="dashed", lwd=2, col="green")
  par(mfrow=c(1,1))
}

### SIRT1
# Additive: None
# Interactive: ChrX, Chr17, Chr14
### BMAL1
# Additive: None
# Interactive: Chr16, Chr13, Chr10, Chr9
### PIK3CG
# Additive: Chr12, ChrX
# Interactive: Chr12, Chr1, Chr9
### NR1H3
# Additive: None
# Interactive: ChrX, Chr16, Chr5, Chr1
### NCOR2
# Additive: None
# Interactive: Chr9, Chr14, Chr19
### SREBF1
# Additive: None
# Interactive: Chr9, Chr16, ChrX

graphics.off()
f2g.scanAdd.scans <- c()
f2g.scanInt.scans <- c()
for (i in 5:length(f2g$pheno)) {
  col <- i - 3
  f2g.scanAdd <- scanone(f2g, pheno.col=4, addcovar=f2g$pheno[,i], method="hk")
  f2g.scanInt <- scanone(f2g, pheno.col=4, intcovar=f2g$pheno[,i], method="hk")
  f2g.scanAdd.scans[names(f2g$pheno)[i]] <- f2g.scanAdd
  f2g.scanInt.scans[names(f2g$pheno)[i]] <- f2g.scanInt
  chrs <- c()
  if (col == 2) {
    chrs <- c("14", "17", "X")
  }
  if (col == 3) {
    chrs <- c("9", "10", "13", "16")
  }
  if (col == 4) {
    chrs <- c("1", "9", "12", "X")
  }
  if (col == 5) {
    chrs <- c("1", "5", "16", "X")
  }
  if (col == 6) {
    chrs <- c("9", "14", "19")
  }
  if (col == 7) {
    chrs <- c("9", "16", "X")
  }
  quartz()
  par(mfrow=c(3,1))
  plot(f2g.scan1, lodcolumn=1, main="FBN1", chr=chrs)
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.63, lty="dashed", lwd=2, col="green")
  plot(f2g.scan1, f2g.scanAdd, lodcolumn=1, main=paste("FBN1 with Additive Covariate", names(f2g$pheno)[i]), chr=chrs)
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.63, lty="dashed", lwd=2, col="green")
  plot(f2g.scan1, f2g.scanInt, lodcolumn=1, main=paste("FBN1 with Interactive Covariate", names(f2g$pheno)[i]), chr=chrs)
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(f2g.scan1, perms=PERMS[,1], alpha=0.63, lty="dashed", lwd=2, col="green")
  par(mfrow=c(1,1))
  print(names(f2g$pheno)[i])
  print(summary(f2g.scanAdd, perms=PERMS, alpha=0.10, format="tabByCol", ci.function="lodint"))
  print(summary(f2g.scanInt, perms=PERMS, alpha=0.10, format="tabByCol", ci.function="lodint"))
}

graphics.off()
### Common FBN1 Chr3 and Chr7 Peaks
quartz()
par(mfrow=c(1,2))
plotPXG(f2g, find.marker(f2g, 3, 49.12), pheno.col="FBN1", main="Chr3 Peak")
plotPXG(f2g, find.marker(f2g, 7, 6.14), pheno.col="FBN1", main="Chr7 Peak")
par(mfrow=c(1,1))
### SIRT1 Interactive ChrX Peak
quartz()
plotPXG(f2g, find.marker(f2g, "X", 7.98), pheno.col="FBN1", main="SIRT1 Interactive ChrX Peak")
### PIK3CG Interactive Chr1 and Chr12 Peaks
quartz()
par(mfrow=c(1,2))
plotPXG(f2g, find.marker(f2g, 1, 78.64), pheno.col="FBN1", main="PIK3CG Interactive Chr1 Peak")
plotPXG(f2g, find.marker(f2g, 12, 21.24), pheno.col="FBN1", main="PIK3CG Interactive Chr12 Peak")
par(mfrow=c(1,1))
### NR1H3 Interactive Chr1, Chr5, and ChrX Peaks
quartz()
par(mfrow=c(1,3)) 
plotPXG(f2g, find.marker(f2g, 1, 74.64), pheno.col="FBN1", main="NR1H3 Interactive Chr1 Peak")
plotPXG(f2g, find.marker(f2g, 5, 3.80), pheno.col="FBN1", main="NR1H3 Interactive Chr5 Peak")
plotPXG(f2g, find.marker(f2g, "X", 28.30), pheno.col="FBN1", main="NR1H3 Interactive ChrX Peak")
par(mfrow=c(1,1))
### NCOR2 Interactive Chr9, Chr14, and Chr19 Peaks
quartz()
par(mfrow=c(1,3))
plotPXG(f2g, find.marker(f2g, 9, 36.00), pheno.col="FBN1", main="NCOR2 Interactive Chr9 Peak")
plotPXG(f2g, find.marker(f2g, 14, 30.92), pheno.col="FBN1", main="NCOR2 Interactive Chr14 Peak")
plotPXG(f2g, find.marker(f2g, 19, 25.56), pheno.col="FBN1", main="NCOR2 Interactive Chr19 Peak")
par(mfrow=c(1,1))
### SREBF1 Interactive Chr1 and Chr9 Peaks
quartz()
par(mfrow=c(1,2))
plotPXG(f2g, find.marker(f2g, 1, 88.35), pheno.col="FBN1", main="SREBF1 Interactive Chr1 Peak")
plotPXG(f2g, find.marker(f2g, 9, 47.29), pheno.col="FBN1", main="SREBF1 Interactive Chr9 Peak")
par(mfrow=c(1,1))

quartz()
par(mfrow=c(1,2))
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 3, 49.12), main="FBN1 Chr3 Common Peak")
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 7, 6.14), main="FBN1 Chr7 Common Peak")
par(mfrow=c(1,1))
quartz()
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, "X", 7.98), main="SIRT1 Interactive ChrX Peak")
quartz()
par(mfrow=c(1,2))
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 1, 78.64), main="PIK3CG Interactive Chr1 Peak")
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 12, 21.24), main="PIK3CG Interactive Chr12 Peak")
par(mfrow=c(1,1))
quartz()
par(mfrow=c(1,3)) 
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 1, 74.64), main="NR1H3 Interactive Chr1 Peak")
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 5, 3.80), main="NR1H3 Interactive Chr5 Peak")
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, "X", 28.30), main="NR1H3 Interactive ChrX Peak")
par(mfrow=c(1,1))
quartz()
par(mfrow=c(1,3))
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 9, 36.00), main="NCOR2 Interactive Chr9 Peak")
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 14, 30.92), main="NCOR2 Interactive Chr14 Peak")
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 19, 25.56), main="NCOR2 Interactive Chr19 Peak")
par(mfrow=c(1,1))
quartz()
par(mfrow=c(1,2))
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 1, 88.35), main="SREBF1 Interactive Chr1 Peak")
effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 9, 47.29), main="SREBF1 Interactive Chr9 Peak")
par(mfrow=c(1,1))

effectplot(f2g, pheno.col="FBN1", mname1=find.marker(f2g, 14, 32.0), mname2=find.marker(f2g, 6, 28.0))
