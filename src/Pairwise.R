#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

source("src/configuration.R")

summary(f2g)

source("src/QTLNet.R")

pairs.simple <- function(data, genes, main="Pair Matrix Scan") {
  
  f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
  for (i in 1:length(genes)) {
    f2g$pheno <- cbind(f2g$pheno, data[,annot$a_gene_id[which(annot$gene_symbol==genes[i])]])
  }
  names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), genes)
  names(f2g$pheno)
  
  # Scatterplot Matrices
  graphics.open()
  pairs(f2g$pheno[,4:length(f2g$pheno)], main=main, labels=genes, upper.panel=panel.cor, diag.panel=panel.hist)
}
pairs.double <- function(data1, genes1, data2, genes2, main="Pair Matrix Scan") {
  
  f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
  for (i in 1:length(genes1)) {
    f2g$pheno <- cbind(f2g$pheno, data1[,annot$a_gene_id[which(annot$gene_symbol==genes1[i])]])
  }
  for (i in 1:length(genes2)) {
    f2g$pheno <- cbind(f2g$pheno, data2[,annot$a_gene_id[which(annot$gene_symbol==genes2[i])]])
  }
  names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), genes1, genes2)
  names(f2g$pheno)
  
  # Scatterplot Matrices
  graphics.open()
  pairs(f2g$pheno[,4:length(f2g$pheno)], main=main, labels=c(genes1, genes2), upper.panel=panel.cor, diag.panel=panel.hist)
}
pairs.phenos.simple <- function(data.g, genes, data.p, phenos, main="Pair Matrix Scan") {
  
  f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
  for (i in 1:length(genes)) {
    f2g$pheno <- cbind(f2g$pheno, data.g[,annot$a_gene_id[which(annot$gene_symbol==genes[i])]])
  }
  for (i in 1:length(phenos)) {
    f2g$pheno <- cbind(f2g$pheno, data.p[,phenos[i]])
  }
  names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), genes, phenos)
  names(f2g$pheno)
  
  # Scatterplot Matrices
  graphics.open()
  pairs(f2g$pheno[,4:length(f2g$pheno)], main=main, upper.panel=panel.cor, diag.panel=panel.hist)
}
pairs.phenos.double <- function(data1.g, genes1, data2.g, genes2, data.p, phenos, main="Pair Matrix Scan") {
  
  f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")] 
  for (i in 1:length(genes1)) {
    f2g$pheno <- cbind(f2g$pheno, data1.g[,annot$a_gene_id[which(annot$gene_symbol==genes1[i])]])
  }
  for (i in 1:length(genes2)) {
    f2g$pheno <- cbind(f2g$pheno, data2.g[,annot$a_gene_id[which(annot$gene_symbol==genes2[i])]])
  }
  for (i in 1:length(phenos)) {
    f2g$pheno <- cbind(f2g$pheno, data.p[,phenos[i]])
  }
  names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), genes1, genes2, phenos)
  names(f2g$pheno)
  
  # Scatterplot Matrices
  graphics.open()
  pairs(f2g$pheno[,4:length(f2g$pheno)], main=main, upper.panel=panel.cor, diag.panel=panel.hist)
}

#################################################
## RECEPTOR STUDIES                            ##
#################################################

# Chromosome 3
genes <- c("Gnai3", "Gnat2", "Gpr61", "Celsr2", "Gpsm2", "Gpr88")
pairs.simple(data=adipose.rz, genes=c(c("Fbn1"), genes), main="FBN1 with Chr3 Adipose Gene Expression")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=liver.rz, genes2=genes, main="FBN1 with Chr3 Liver Gene Expression")
pairs.simple(data=liver.rz, genes=c(c("Fbn1"), genes), main="FBN1 (Liver) with Chr3 Liver Gene Expression")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=adipose.rz, genes2=genes, main="FBN1 (Liver) with Chr3 Adipose Gene Expression")

# Chromosome 7
genes <- c("Grik5", "Gng8", "Gpr4")
pairs.simple(data=adipose.rz, genes=c(c("Fbn1"), genes), main="FBN1 with Chr7 Adipose Gene Expression")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=liver.rz, genes2=genes, main="FBN1 with Chr7 Liver Gene Expression")
pairs.simple(data=liver.rz, genes=c(c("Fbn1"), genes), main="FBN1 (Liver) with Chr7 Liver Gene Expression")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=adipose.rz, genes2=genes, main="FBN1 (Liver) with Chr7 Adipose Gene Expression")

# Chromosome 8
genes <- c("Npy1r", "Npy5r", "Tm6sf2", "Gpr114")
pairs.simple(data=adipose.rz, genes=c(c("Fbn1"), genes), main="FBN1 with Chr8 Adipose Gene Expression")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=liver.rz, genes2=genes, main="FBN1 with Chr8 Liver Gene Expression")
pairs.simple(data=liver.rz, genes=c(c("Fbn1"), genes), main="FBN1 (Liver) with Chr8 Liver Gene Expression")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=adipose.rz, genes2=genes, main="FBN1 (Liver) with Chr8 Adipose Gene Expression")

# GNG8-Related Pathway: GHRH AKT Signaling
genes <- c("Akt3", "Pfkfb2", "Hk2", "Gys1", "Gys2", "Tbc1d4", "Acly", "Cd19")
pairs.simple(data=adipose.rz, genes=c(c("Fbn1"), genes), main="Adipose Gene Expression")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=liver.rz, genes2=genes, main="Liver Gene Expression")
pairs.simple(data=liver.rz, genes=c(c("Fbn1"), genes), main="Liver Gene Expression (w/ FBN1 Liver)")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=adipose.rz, genes2=genes, main="Adipose Gene Expression (w/ FBN1 Liver)")

# GNG8-Related Pathway: GHRH PKA Signaling
genes <- c("Prkaca", "Prkacb", "Pygm", "Lipe")
pairs.simple(data=adipose.rz, genes=c(c("Fbn1"), genes), main="Adipose Gene Expression")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=liver.rz, genes2=genes, main="Liver Gene Expression")
pairs.simple(data=liver.rz, genes=c(c("Fbn1"), genes), main="Liver Gene Expression (w/ FBN1 Liver)")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=adipose.rz, genes2=genes, main="Adipose Gene Expression (w/ FBN1 Liver)")

# PKA/PKB/GPCR Signaling
genes <- c("Gpr12", "Rabgap1,Gpr21", "Gpr39", "Gpr50")
genes <- c("Gpr61", "Gpr62", "Gpr82", "Gpr119")
genes <- c("Gpr146", "Gpr171", "Gprc5b", "Tpra1")
pairs.simple(data=adipose.rz, genes=c(c("Fbn1"), genes), main="Adipose Gene Expression")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=liver.rz, genes2=genes, main="Liver Gene Expression")
pairs.simple(data=liver.rz, genes=c(c("Fbn1"), genes), main="Liver Gene Expression (w/ FBN1 Liver)")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=adipose.rz, genes2=genes, main="Adipose Gene Expression (w/ FBN1 Liver)")

# Class A Orphan Receptor Candidates
# Gpr12, Gpr21, Gpr50, Gpr62, Gpr82, Gpr119, Gpr146, Gprc5b
genes <- c("Gpr12", "Rabgap1,Gpr21", "Gpr50", "Gpr62", "Gpr82", "Gpr119", "Gpr146", "Gprc5b")
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("INS.4wk", "NEFA")]), phenos=c("INS.4wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("INS.6wk", "NEFA")]), phenos=c("INS.6wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("INS.8wk", "NEFA")]), phenos=c("INS.8wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("INS.10wk", "NEFA")]), phenos=c("INS.10wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("GLU.4wk", "NEFA")]), phenos=c("GLU.4wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("GLU.6wk", "NEFA")]), phenos=c("GLU.6wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("GLU.8wk", "NEFA")]), phenos=c("GLU.8wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("GLU.10wk", "NEFA")]), phenos=c("GLU.10wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("HOMA.4wk", "NEFA")]), phenos=c("HOMA.4wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("HOMA.6wk", "NEFA")]), phenos=c("HOMA.6wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("HOMA.8wk", "NEFA")]), phenos=c("HOMA.8wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("HOMA.10wk", "NEFA")]), phenos=c("HOMA.10wk"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("Liver.wt.1", "NEFA")]), phenos=c("Liver.wt.1"))
pairs.phenos.simple(data.g=liver.rz, genes=genes, data.p=log(phenotypes[,c("Fat.wt", "NEFA")]), phenos=c("Fat.wt"))

# Class A Orphan Receptor Candidates
# Gpr12, Gpr21, Gprc5b

runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr12", "Prkaca", "Akt3"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr12"), phenos=c("INS.10wk"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr12"), phenos=c("GLU.10wk"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr12"), phenos=c("Fat.wt"))

annot["10028674547.1",]$gene_symbol <- "Gpr21"
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr21", "Prkaca", "Akt3"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr21"), phenos=c("INS.10wk"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr21"), phenos=c("GLU.10wk"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gpr21"), phenos=c("Fat.wt"))
annot["10028674547.1",]$gene_symbol <- "Rabgap1,Gpr21"

runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gprc5b", "Prkaca", "Akt3"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gprc5b"), phenos=c("INS.10wk"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gprc5b"), phenos=c("GLU.10wk"))
runQTLNet(genes.adipose=c("Fbn1"), genes.liver=c("Gprc5b"), phenos=c("Fat.wt"))

# Insulin-Transduction-Like Tyrosine Kinase Signaling
genes <- c("Insr", "Igf1r")
genes <- c("Fgfr1", "Fgfr3", "Fgfr4")
genes <- c("Pik3ca", "Pik3cb", "Pik3cg", "Gsk3b")
pairs.simple(data=adipose.rz, genes=c(c("Fbn1"), genes), main="Adipose Gene Expression")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=liver.rz, genes2=genes, main="Liver Gene Expression")
pairs.simple(data=liver.rz, genes=c(c("Fbn1"), genes), main="Liver Gene Expression (w/ FBN1 Liver)")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=adipose.rz, genes2=genes, main="Adipose Gene Expression (w/ FBN1 Liver)")

# Previous research has shown that Asprosin acts as an effector adipokine
# Tests have revealed that it binds to surface receptors on hepatocytes
# Further tests have shown that it utilizes the G Protein-cAMP-PKA axis to do so
# Insulin was also found to suppress the action of asprosin (competitively)

# Our initial tests suggest that Gng8 (G Protein gamma, part of the beta-gamma machinery) is correlated with Fbn1 expression
# Since Insulin competitevly inhibits asprosin action, there is also the likelihood of a tyrosine-kinase receptor
# We identify Class A G-Protein Coupled Receptors (Orphans) that are correlated with Fbn1
# Results also suggest that the pathways either overlap or a tyrosine-kinase receptor may also affect the situation

# Pancreatic / Islet Tissue
genes <- c("Npy1r", "Npy5r", "Tm6sf2", "Gpr114")
pairs.double(data1=adipose.rz, genes1=c("Fbn1"), data2=islet.rz, genes2=genes, main="Islet Gene Expression")
pairs.phenos.simple(data.g=islet.rz, genes=genes, data.p=phenotypes.rz, phenos=c("INS.10wk", "GLU.4wk", "HOMA.8wk"))

# Chromosome 2
genes <- c("Pcsk2", "Rrbp1", "Ralgapa2", "Insm1", "Foxa2", "Acss1", "Trib3")
pairs.double(data1=liver.rz, genes1=c("Fbn1"), data2=islet.rz, genes2=genes, main="Islet Gene Expression")
pairs.phenos.simple(data.g=islet.rz, genes=genes, data.p=phenotypes.rz, phenos=c("INS.10wk", "GLU.4wk", "HOMA.8wk"))
pairs.double(data1=liver.rz, genes1=c("Fbn1", "Gprc5b", "Rabgap1,Gpr21", "Trib3"), data2=islet.rz, genes2=c("Trib3"))

pairs.phenos.simple(data.g=adipose.rz, genes=c("Amy1", "Agl", "Fabp2", "Cyp2u1", "Fbn1"), data.p=phenotypes.rz, phenos=c("HOMA.10wk"))
pairs.phenos.simple(data.g=adipose.rz, genes=c("Lipe", "Cyp2s1", "Cyp2b10", "Foxa3", "Apoc4", "Apoc2", "Apoe", "Gsk3a", "Fbn1"), data.p=phenotypes.rz, phenos=c("HOMA.10wk"))
pairs.phenos.simple(data.g=adipose.rz, genes=c("Amy1", "Cyp2s1", "Apoc4", "Apoc2", "Apoe", "Fbn1"), data.p=phenotypes.rz, phenos=c("HOMA.10wk"))
