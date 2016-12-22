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
  x11()
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
  x11()
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
  x11()
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
  x11()
  pairs(f2g$pheno[,4:length(f2g$pheno)], main=main, upper.panel=panel.cor, diag.panel=panel.hist)
}