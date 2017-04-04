#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: March 20 2016                         ##
#################################################

source("src/configuration.R")
source("src/ScanOne.R")

#################################################
## FULL SCRUB                                  ##
#################################################

THRESHOLD <- 0.50
INS.THRESHOLD <- 0.20
HOMA.THRESHOLD <- 0.20
FAT.THRESHOLD <- 0.20

genes <- c()
corrs <- c()
for (i in 1:nrow(annot)) {
  
  gene <- annot[i,"gene_symbol"]
  if (gene.data.exists(gene)) {
    data1 <- gene.exp("Sirt2", adipose.rz)
    data.ins <- clinical("INS.10wk", phenotypes.rz)
    data.homa <- clinical("HOMA.8wk", phenotypes.rz)
    data.fat <- clinical("Fat.wt", phenotypes.rz)
    data2 <- gene.exp(gene, adipose.rz)
    indx <- sort(unique(c(which(is.na(data1)), which(is.na(data.ins)), which(is.na(data.homa)), which(is.na(data.fat)), which(is.na(data2)))))
    data1 <- data1[-indx]
    data.ins <- data.ins[-indx]
    data.homa <- data.homa[-indx]
    data.fat <- data.fat[-indx]
    data2 <- data2[-indx]
    correlation <- abs(cor(x=data1, y=data2))
    correlation.ins <- abs(cor(x=data.ins, y=data2))
    correlation.homa <- abs(cor(x=data.homa, y=data2))
    correlation.fat <- abs(cor(x=data.fat, y=data2))
    if ((correlation >= THRESHOLD) && (correlation.ins >= INS.THRESHOLD && correlation.homa >= HOMA.THRESHOLD && correlation.fat >= FAT.THRESHOLD)) {
      print(paste("Gene: ", gene, " | Correlation: ", correlation, " | Insulin: ", correlation.ins, " | HOMA: ", correlation.homa, " | Fat.wt: ", correlation.fat, sep=""))
      genes <- c(genes, paste(gene))
      corrs <- c(corrs, correlation)
    }
  }
}

gene.data.frame <- data.frame(genes, corrs)
write.csv(file="test.csv", gene.data.frame)

#################################################
## PARTIAL SCRUB                               ##
#################################################

CSV <- "data/chr-7-genes.csv"
PDF <- "report/paper/genes/Ins-0.30-ins-0.20-homa-0.20-fat-0.20.pdf"
THRESHOLD <- 0.20

chrs <- read.csv(CSV)
genes <- c()
for (i in 1:nrow(chrs)) {
  
  gene <- chrs[i,1]
  if (gene.data.exists(gene)) {
    data1 <- gene.exp("Fbn1", adipose.rz)
    data2 <- gene.exp(gene, adipose.rz)
    indx <- sort(unique(c(which(is.na(data1)), which(is.na(data2)))))
    data1 <- data1[-indx]
    data2 <- data2[-indx]
    correlation <- cor(x=data1, y=data2)
    if (correlation >= THRESHOLD) {
      print(paste("Gene: ", gene, " | Correlation: ", correlation, sep=""))
      genes <- c(genes, paste(gene))
    }
  }  
}

run <- function() {
  pdf(PDF)
  graphics.allow.open <- FALSE
  scan.fbn1 <- runScanOne(clinical("INS.10wk", phenotypes.rz), "INS.10wk.perms.RData", perms=100, clean=TRUE, graphs=FALSE)
  for (i in 1:length(genes)) {
    
    print(genes[i])
    name <- genes[i]
    covar <- gene.exp(genes[i], adipose.rz)
    
    scan.fbn1.add.apoe <- runScanOne.AddCovar(clinical("INS.10wk", phenotypes.rz), "INS.10wk.perms.RData", covar, 100, clean=TRUE, graphs=FALSE)
    scan.fbn1.int.apoe <- runScanOne.IntCovar(clinical("INS.10wk", phenotypes.rz), "INS.10wk.perms.RData", covar, 100, clean=TRUE, graphs=FALSE)
    
    par(mfrow=c(3,1))
    
    plot(scan.fbn1$scan, lodcolumn=1, main="INS.10wk LOD for Variation in Mouse Genome")
    add.threshold(scan.fbn1$scan, perms=scan.fbn1$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
    add.threshold(scan.fbn1$scan, perms=scan.fbn1$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
    add.threshold(scan.fbn1$scan, perms=scan.fbn1$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
    
    plot(scan.fbn1.add.apoe$scan, lodcolumn=1, main=paste("INS.10wk with ", name, " as Additive Covariate", sep=""))
    add.threshold(scan.fbn1.add.apoe$scan, perms=scan.fbn1.add.apoe$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
    add.threshold(scan.fbn1.add.apoe$scan, perms=scan.fbn1.add.apoe$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
    add.threshold(scan.fbn1.add.apoe$scan, perms=scan.fbn1.add.apoe$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
    
    plot(scan.fbn1.int.apoe$scan, lodcolumn=1, main=paste("INS.10wk with ", name, " as Interactive Covariate", sep=""))
    add.threshold(scan.fbn1.int.apoe$scan, perms=scan.fbn1.int.apoe$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
    add.threshold(scan.fbn1.int.apoe$scan, perms=scan.fbn1.int.apoe$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
    add.threshold(scan.fbn1.int.apoe$scan, perms=scan.fbn1.int.apoe$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
    
    par(mfrow=c(1,1))
  }
  dev.off()
  graphics.allow.open <- TRUE
}

run()