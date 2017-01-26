#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: January 12 2016                       ##
#################################################

source("src/configuration.R")
source("src/ScanOne.R")

# symbol: The gene symbol to be checked
gene.data.exists <- function(symbol) {
  
  return(length(annot$gene_symbol[which(annot$gene_symbol==symbol)]) > 0)
}

# gene.name: The common MGI gene symbol
# data.set: The gene expression dataset, such as adipose.rz
gene.exp <- function(gene.name, data.set) {
  
  if (gene.data.exists(gene.name)) {
    
    return(data.set[,annot$a_gene_id[which(annot$gene_symbol==gene.name)]])
  }
}

# clin.name: The name of the clinical trait in the dataset
# data.set: Either phenotypes or phenotypes.rz
clinical <- function(clin.name, data.set) {
  
  if (clin.name %in% names(data.set)) {
    
    return(data.set[,clin.name])
  }
}

# chr: The chromosome on which to look for the genotype data
# pos: The position at which to look for genotype data
genotype <- function(chr, pos) {
  
  if (chr > 0 && chr < 21 && pos >= 0 && pos <= 100) {
    
    return(f2g$geno[[chr]]$data[,find.marker(f2g, chr=chr, pos=pos)])
  }
}

# Run a scanone on Fbn1
Fbn1 <- gene.exp(gene.name="Fbn1", data.set=adipose.rz)
runScanOne(params=Fbn1, file.name="FBN1.adipose.perms.RData", perms=100)

#[1] "MouseNum" "Sex"      "pgm"      "params"  
#           chr   pos ci.low ci.high  lod
#rs13477374   3 49.12  44.55    65.1 4.40
#c7.loc4      7  6.14   2.14    15.9 4.07
#           chr   pos ci.low ci.high  lod
#c1.loc75     1 76.64   1.64    86.6 2.83
#rs13477374   3 49.12  44.55    65.1 4.40
#c7.loc4      7  6.14   2.14    15.9 4.07
#cX.loc5      X  8.60   3.60    28.3 2.17
#[1] "Chromosome  3"
#[1] 46.54620 63.65573
#[1] "Chromosome  7"
#[1]  2.1403 15.1403

# The QTL on Chromosome 3, the effector
QTL3 <- genotype(chr=3, pos=49.12)
# Glucose is affected
GLU.8wk <- log(clinical(clin.name="GLU.8wk", data.set=phenotypes))
# Fbn1 (Asprosin) is the mediator 
Fbn1 <- gene.exp(gene.name="Fbn1", data.set=adipose.rz)
# Sex is the covariate
Sex <- as.numeric(f2g$pheno[,"Sex"]) - 1

# Remove all NA values from the data
indx <- sort(unique(c(which(is.na(QTL3)), which(is.na(GLU.8wk)), which(is.na(Fbn1)))))
QTL3 <- QTL3[-indx]
GLU.8wk <- GLU.8wk[-indx]
Fbn1 <- Fbn1[-indx]
Sex <- Sex[-indx]

# Create a data frame to house the data we will work with
model.data <- data.frame(QTL3, GLU.8wk, Fbn1, Sex)

# This is the fit model for the mediator
med.fit <- lm(Fbn1 ~ Sex + QTL3, data=model.data)
# This is the fit model for the affected
out.fit <- glm(GLU.8wk ~ Sex + Fbn1 + QTL3, data=model.data)

# The output of the mediation analysis
med.out <- mediate(med.fit, out.fit, treat="QTL3", mediator="Fbn1", robustSE=TRUE, sims=100)

# Print summary and view scores
summary(med.out)
plot(med.out)

