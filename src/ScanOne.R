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