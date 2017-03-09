#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

source("src/configuration.R")
source("src/ScanOne.R")

CSV <- "data/chr-7-genes.csv"
PDF <- "report/paper/genes/chr3.pdf"
THRESHOLD <- 0.15

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

pdf(PDF)
graphics.allow.open <- FALSE
scan.fbn1 <- runScanOne(gene.exp("Fbn1", adipose.rz), "Fbn1.adipose.perms.RData", perms=100, clean=TRUE, graphs=FALSE)
for (i in 1:length(genes)) {
  
  name <- genes[i]
  covar <- gene.exp(genes[i], adipose.rz)
  
  scan.fbn1.add.apoe <- runScanOne.AddCovar(gene.exp("Fbn1", adipose.rz), "Fbn1.adipose.perms.RData", covar, 100, clean=TRUE, graphs=FALSE)
  scan.fbn1.int.apoe <- runScanOne.IntCovar(gene.exp("Fbn1", adipose.rz), "Fbn1.adipose.perms.RData", covar, 100, clean=TRUE, graphs=FALSE)
  
  par(mfrow=c(3,1))
  
  plot(scan.fbn1$scan, lodcolumn=1, main="FBN1 LOD for Variation in Mouse Genome")
  add.threshold(scan.fbn1$scan, perms=scan.fbn1$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
  add.threshold(scan.fbn1$scan, perms=scan.fbn1$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
  add.threshold(scan.fbn1$scan, perms=scan.fbn1$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
  
  plot(scan.fbn1.add.apoe$scan, lodcolumn=1, main=paste("FBN1 with ", name, " as Additive Covariate", sep=""))
  add.threshold(scan.fbn1.add.apoe$scan, perms=scan.fbn1.add.apoe$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
  add.threshold(scan.fbn1.add.apoe$scan, perms=scan.fbn1.add.apoe$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
  add.threshold(scan.fbn1.add.apoe$scan, perms=scan.fbn1.add.apoe$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
  
  plot(scan.fbn1.int.apoe$scan, lodcolumn=1, main=paste("FBN1 with ", name, " as Interactive Covariate", sep=""))
  add.threshold(scan.fbn1.int.apoe$scan, perms=scan.fbn1.int.apoe$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
  add.threshold(scan.fbn1.int.apoe$scan, perms=scan.fbn1.int.apoe$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
  add.threshold(scan.fbn1.int.apoe$scan, perms=scan.fbn1.int.apoe$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
  
  par(mfrow=c(1,1))
}
dev.off()
graphics.allow.open <- TRUE