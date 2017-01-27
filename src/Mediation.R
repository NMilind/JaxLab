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

Fbn1 <- gene.exp("Fbn1", data.set=adipose.rz)
GLU.4wk <- log(clinical(clin.name="GLU.4wk", data.set=phenotypes))
Gpr21 <- gene.exp(gene.name="Rabgap1,Gpr21", data.set=liver.rz)

# Remove all NA values from the data
indx <- sort(unique(c(which(is.na()))))

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
