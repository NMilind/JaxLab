#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

#################################################
## 1: ENVIRONMENTAL SETUP                      ##
#################################################

source("src/configuration.R")

summary(f2g)

# Proteins Tested Before QTL Analysis
#genes.adipose <- c("Fbn1", "Sirt1", "Pik3cg", "Nr1h3", "Ncor2", "Srebf1", "Tgfb1")
# Proteins on Chr3 QTL Confidence Interval
genes.adipose <- c("Fbn1", "Gpr61", "Gpsm2", "Fndc7", "Col11a1", "S1pr1", "Gpr88", "Agl", "Cnn3", "Myoz2", "Pla2g12a", "Col25a1", "Cyp2u1", "Emcn")
genes.liver <- c("Fbn1", "Gpr61", "Gpsm2", "Fndc7", "Col11a1", "S1pr1", "Gpr88", "Agl", "Cnn3", "Myoz2", "Pla2g12a", "Col25a1", "Cyp2u1", "Emcn")
# Proteins on Chr7 QTL Confidence Interval
#genes.adipose <- c("Fbn1", "Cdc42ep5", "Shisa7", "Gng8", "Gpr4", "Kptn", "Mark4", "Strn4", "Tgfb1")
#genes.liver <- c("Fbn1", "Cdc42ep5", "Shisa7", "Gng8", "Gpr4", "Kptn", "Mark4", "Strn4", "Tgfb1")

f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
for (i in 1:length(genes.adipose)) {
  f2g$pheno <- cbind(f2g$pheno, adipose.rz[,annot$a_gene_id[which(annot$gene_symbol==genes.adipose[i])]])
}
names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), genes.adipose)
names(f2g$pheno)

# Scatterplot Matrices
x11()
pairs(f2g$pheno[,4:length(f2g$pheno)], main="Adipose", labels=genes.adipose, upper.panel=panel.cor, diag.panel=panel.hist)

f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
for (i in 1:length(genes.liver)) {
  f2g$pheno <- cbind(f2g$pheno, liver.rz[,annot$a_gene_id[which(annot$gene_symbol==genes.liver[i])]])
}
names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), genes.liver)
names(f2g$pheno)

# Scatterplot Matrices
x11()
pairs(f2g$pheno[,4:length(f2g$pheno)], main="Liver", labels=genes.liver, upper.panel=panel.cor, diag.panel=panel.hist)

# Interesting Chr3 QTL Genes
#   Gpsm2: 0.40
#   Cnn3: 0.30
#   Pla2g12a: -0.26
# Interesting Chr7 QTL Genes
#   Gpsm2: 0.30
#   S1pr1: 0.33
#   Agl: 0.23
#   Pla2g12a: -0.37
#   Cyp2u1: -0.46
#   Emcn: 0.37