#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

source("src/configuration.R")

summary(f2g)

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
  pairs(f2g$pheno[,4:length(f2g$pheno)], main=main, labels=cbind(genes1, genes2), upper.panel=panel.cor, diag.panel=panel.hist)
}

