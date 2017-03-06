#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

runScanOne <- function(params, file.name, perms=1000) {
  
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
  f2g$pheno <- cbind(f2g$pheno[,names(f2g$pheno)], params)
  #  Check to make sure param is added
  print(names(f2g$pheno))
  
  #################################################
  ## SCANS                                       ##
  #################################################
  
  f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
  f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
  
  sex <- as.numeric(f2g$pheno$Sex) - 1
  
  # Permute through the cross to create baseline LOD score levels
  if (!file.exists(paste("data/", file.name, sep=""))) {
    perms <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk", n.perm=perms, perm.Xsp=TRUE)
    save(file=paste("data/", file.name, sep=""), perms)
  }
  else {
    load(file=paste("data/", file.name, sep=""))
  }
  
  # Scan for param in the cross
  scan1 <- scanone(f2g, pheno.col=4, addcovar=sex, method="hk")
  
  # Plot the LOD scan with thresholds from the permutations
  graphics.open()
  plot(scan1, lodcolumn=1, main="LOD for Variation across Mouse Genome")
  add.threshold(scan1, perms=perms, alpha=0.05, lty="dashed", lwd=1, col="green")
  add.threshold(scan1, perms=perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
  add.threshold(scan1, perms=perms, alpha=0.63, lty="dashed", lwd=1, col="red")
  
  # Tabulate and print LOD peaks where alpha=0.05
  print(summary(scan1, perms=perms, alpha=0.05, format="tabByCol", ci.function="lodint"))
  major.peaks <- summary(scan1, perms=perms, alpha=0.05)
  
  # Tabulate and print LOD peaks where alpha=0.63
  print(summary(scan1, perms=perms, alpha=0.63, format="tabByCol", ci.function="lodint"))
  
  if (nrow(major.peaks) != 0) {
    for (i in seq(nrow(major.peaks))) {
      graphics.open()
      chr <- major.peaks[i,1]
      ci <- bayesint(scan1, chr=chr, prob=0.95)
      plot(scan1, chr=chr, lodcolumn=1, main=paste("Confidence Interval Chr", chr, sep=""))
      lines(x=ci[c(1,3),2], y=c(0,0), type="l", col="#00FF00", lwd=4)
      print(paste("Chromosome ", chr))
      print(ci[c(1,3),2])
    }
  }
  
  return(list("scan"=scan1, "perms"=perms))
}

#scan.amy1 <- runScanOne(gene.exp("Amy1", data=adipose.rz), file.name="Amy1.adipose.perms.RData", perms=100)
#scan.cyp2s1 <- runScanOne(gene.exp("Cyp2s1", data=adipose.rz), file.name="Cyp2s1.adipose.perms.RData", perms=100)
#scan.apoc4 <- runScanOne(gene.exp("Apoc4", data=adipose.rz), file.name="Apoc4.adipose.perms.RData", perms=100)
#scan.apoc2 <- runScanOne(gene.exp("Apoc2", data=adipose.rz), file.name="Apoc2.adipose.perms.RData", perms=100)
#scan.apoe <- runScanOne(gene.exp("Apoe", data=adipose.rz), file.name="Apoe.adipose.perms.RData", perms=100)
#scan.homa <- runScanOne(clinical(clin.name="HOMA.10wk", data.set=phenotypes.rz), file.name="HOMA.10wk.perms.RData", perms=100)

#graphics.open()
#par(mfrow=c(3, 2))
#plot(scan.amy1$scan, lodcolumn=1, main="ScanOne Plot for Amy1 in Adipose Tissue")
# add.threshold(scan.amy1$scan, perms=scan.amy1$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
# add.threshold(scan.amy1$scan, perms=scan.amy1$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
# add.threshold(scan.amy1$scan, perms=scan.amy1$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# 
# plot(scan.cyp2s1$scan, lodcolumn=1, main="ScanOne Plot for Cyp2s1 in Adipose Tissue")
# add.threshold(scan.cyp2s1$scan, perms=scan.cyp2s1$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
# add.threshold(scan.cyp2s1$scan, perms=scan.cyp2s1$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
# add.threshold(scan.cyp2s1$scan, perms=scan.cyp2s1$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# 
# plot(scan.apoc4$scan, lodcolumn=1, main="ScanOne Plot for Apoc4 in Adipose Tissue")
# add.threshold(scan.apoc4$scan, perms=scan.apoc4$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
# add.threshold(scan.apoc4$scan, perms=scan.apoc4$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
# add.threshold(scan.apoc4$scan, perms=scan.apoc4$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# 
# plot(scan.apoc2$scan, lodcolumn=1, main="ScanOne Plot for Apoc2 in Adipose Tissue")
# add.threshold(scan.apoc2$scan, perms=scan.apoc2$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
# add.threshold(scan.apoc2$scan, perms=scan.apoc2$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
# add.threshold(scan.apoc2$scan, perms=scan.apoc2$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# 
# plot(scan.apoe$scan, lodcolumn=1, main="ScanOne Plot for Apoe in Adipose Tissue")
# add.threshold(scan.apoe$scan, perms=scan.apoe$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
# add.threshold(scan.apoe$scan, perms=scan.apoe$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
# add.threshold(scan.apoe$scan, perms=scan.apoe$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# 
# plot(scan.homa$scan, lodcolumn=1, main="ScanOne Plot for HOMA(10wk) in Adipose Tissue")
# add.threshold(scan.homa$scan, perms=scan.homa$perms, alpha=0.05, lty="dashed", lwd=1, col="green")
# add.threshold(scan.homa$scan, perms=scan.homa$perms, alpha=0.10, lty="dashed", lwd=1, col="orange")
# add.threshold(scan.homa$scan, perms=scan.homa$perms, alpha=0.63, lty="dashed", lwd=1, col="red")
# 
# par(mfrow=c(1,1))
