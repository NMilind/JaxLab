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

#Fbn1 <- gene.exp(gene.name="Fbn1", data.set=adipose.rz)
Fbn1 <- genotype(chr=3, pos=49.5)
GLU.4wk <- log(clinical(clin.name="GLU.4wk", data.set=phenotypes))
INS.10wk <- log(clinical(clin.name="INS.10wk", data.set=phenotypes))
HOMA.10wk <- (clinical(clin.name="HOMA.10wk", data.set=phenotypes.rz))
Gpr21 <- gene.exp(gene.name="Rabgap1,Gpr21", data.set=liver.rz)
Gprc5b <- gene.exp(gene.name="Gprc5b", data.set=liver.rz)
Akt3 <- gene.exp(gene.name="Akt3", data.set=liver.rz)
Prkaca <- gene.exp(gene.name="Prkaca", data.set=liver.rz)
Tgfb2 <- gene.exp(gene.name="Tgfb2", data.set=adipose.rz)
Timp2 <- gene.exp(gene.name="Timp2", data.set=adipose.rz)
Sex <- as.numeric(f2g$pheno[,"Sex"]) - 1

# Remove all NA values from the data
indx <- sort(unique(c(
  which(is.na(Fbn1)), 
  which(is.na(GLU.4wk)), 
  which(is.na(INS.10wk)), 
  which(is.na(HOMA.10wk)),
  which(is.na(Gpr21)), 
  which(is.na(Gprc5b)), 
  which(is.na(Prkaca)), 
  which(is.na(Akt3))))
)
Fbn1 <- Fbn1[-indx]
GLU.4wk <- GLU.4wk[-indx]
INS.10wk <- INS.10wk[-indx]
HOMA.10wk <- HOMA.10wk[-indx]
Gpr21 <- Gpr21[-indx]
Gprc5b <- Gprc5b[-indx]
Akt3 <- Akt3[-indx]
Prkaca <- Prkaca[-indx]
Sex <- Sex[-indx]

# Create a data frame to house the data we will work with
model.data <- data.frame(Fbn1, GLU.4wk, INS.10wk, HOMA.10wk, Gpr21, Gprc5b, Akt3, Prkaca, Sex)

affected <- HOMA.10wk

par(mfrow=c(2,2))

#################################################
## Gpr21 with PKA                              ##
#################################################

# This is the fit model for the mediator
med.fit <- lm(Prkaca ~ Sex + Gpr21 + Fbn1, data=model.data)
# This is the fit model for the affected
out.fit <- glm(affected ~ Sex + Gpr21 + Prkaca + Fbn1, data=model.data)

# The output of the mediation analysis
med.out <- mediate(med.fit, out.fit, treat="Gpr21", mediator="Prkaca", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Gpr21 with PKA")

#################################################
## Gpr21 with PKB                              ##
#################################################

# This is the fit model for the mediator
med.fit <- lm(Akt3 ~ Sex + Gpr21 + Fbn1, data=model.data)
# This is the fit model for the affected
out.fit <- glm(affected ~ Sex + Gpr21 + Akt3 + Fbn1, data=model.data)

# The output of the mediation analysis
med.out <- mediate(med.fit, out.fit, treat="Gpr21", mediator="Akt3", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Gpr21 with PKB")

#################################################
## Gprc5b with PKA                             ##
#################################################

# This is the fit model for the mediator
med.fit <- lm(Prkaca ~ Sex + Gprc5b + Fbn1, data=model.data)
# This is the fit model for the affected
out.fit <- glm(affected ~ Sex + Gprc5b + Prkaca + Fbn1, data=model.data)

# The output of the mediation analysis
med.out <- mediate(med.fit, out.fit, treat="Gprc5b", mediator="Prkaca", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Gprc5b with PKA")

#################################################
## Gprc5b with PKB                             ##
#################################################

# This is the fit model for the mediator
med.fit <- lm(Akt3 ~ Sex + Gprc5b + Fbn1, data=model.data)
# This is the fit model for the affected
out.fit <- glm(affected ~ Sex + Gprc5b + Akt3 + Fbn1, data=model.data)

# The output of the mediation analysis
med.out <- mediate(med.fit, out.fit, treat="Gprc5b", mediator="Akt3", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Gprc5b with PKB")

par(mfrow=c(1,1))

#################################################
## TIMP2 and NFATC4                            ##
#################################################

#Fbn1 <- gene.exp(gene.name="Fbn1", data.set=adipose.rz)
Fbn1 <- genotype(chr=3, pos=49.5)
INS.10wk <- log(clinical(clin.name="INS.10wk", data.set=phenotypes))
Tgfb2 <- gene.exp(gene.name="Tgfb2", data.set=adipose.rz)
Timp2 <- gene.exp(gene.name="Timp2", data.set=adipose.rz)
Nfatc4 <- gene.exp(gene.name="Nfatc4", data.set=adipose.rz)

# Remove all NA values from the data
indx <- sort(unique(c(
  which(is.na(Fbn1)),
  which(is.na(INS.10wk)),
  which(is.na(Tgfb2)),
  which(is.na(Timp2)),
  which(is.na(Nfatc4))
)))

Fbn1 <- Fbn1[-indx]
INS.10wk <- INS.10wk[-indx]
Tgfb2 <- Tgfb2[-indx]
Timp2 <- Timp2[-indx]
Nfatc4 <- Nfatc4[-indx]

model.data <- data.frame(Fbn1, INS.10wk, Tgfb2, Timp2, Nfatc4)

par(mfrow=c(2,2))

med.fit <- lm(Timp2 ~ Tgfb2 + Fbn1, data=model.data)
out.fit <- glm(INS.10wk ~ Timp2 + Tgfb2 + Fbn1, data=model.data)

med.out <- mediate(med.fit, out.fit, treat="Tgfb2", mediator="Timp2", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Effect of Tgfb2 on Insulin (10 weeks) mediated by Timp2")

med.fit <- lm(Nfatc4 ~ Tgfb2 + Fbn1, data=model.data)
out.fit <- glm(INS.10wk ~ Nfatc4 + Tgfb2 + Fbn1, data=model.data)

med.out <- mediate(med.fit, out.fit, treat="Tgfb2", mediator="Nfatc4", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Effect of Tgfb2 on Insulin (10 weeks) mediated by Nfatc4")

med.fit <- lm(Timp2 ~ Tgfb2 + Fbn1, data=model.data)
out.fit <- glm(INS.10wk ~ Timp2 + Tgfb2 + Fbn1, data=model.data)

med.out <- mediate(med.fit, out.fit, treat="Fbn1", mediator="Timp2", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Effect of Fbn1 on Insulin (10 weeks) mediated by Timp2")

med.fit <- lm(Nfatc4 ~ Tgfb2 + Fbn1, data=model.data)
out.fit <- glm(INS.10wk ~ Nfatc4 + Tgfb2 + Fbn1, data=model.data)

med.out <- mediate(med.fit, out.fit, treat="Fbn1", mediator="Nfatc4", robustSE=TRUE, sims=100)

summary(med.out)
plot(med.out, main="Effect of Fbn1 on Insulin (10 weeks) mediated by Nfatc4")

par(mfrow=c(1,1))
