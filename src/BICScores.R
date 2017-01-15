#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: January 12 2016                       ##
#################################################

source("src/configuration.R")

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

# X: gene expression data for a given gene
# Y: quantitative measurement of clinical phenotype
# Q: genotype at a given marker
triple.fit <- function(X, Y, Q) {
  
  # Remove any NA values from the data
  indx <- sort(unique(c(which(is.na(X)), which(is.na(Y)), which(is.na(Q)))))
  X <- X[-indx]
  Y <- Y[-indx]
  Q <- Q[-indx]
  print(paste("Removed ", length(indx), " rows with NA values from data.", sep=""))
  
  # Calculate BIC scores for models
  bic.independent <- BIC(lm(X~Q)) + BIC(lm(Y~Q)) # X<-Q->Y
  bic.reactive <- BIC(lm(X~Y)) + BIC(lm(Y~Q)) # Q->Y->X
  bic.causal <- BIC(lm(X~Q)) + BIC(lm(Y~X)) # Q->X->Y
  bic.complex <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X))
  
  # Print out the scores from each model
  #print("BIC Scores of each model")
  scores <- c(bic.independent, bic.reactive, bic.causal, bic.complex)
  names(scores) <- c("independent", "reactive", "causal", "complex")
  #print(scores)
  
  # Make lowest BIC score 0 and linearize all other scores accordingly to calculate Delta values
  deltas <- scores - min(scores)
  
  # Print delta values
  #print("Delta values of BIC Scores - should be larger than 10 for significance")
  #print(deltas)
  
  # Estimate the strength of evidence for each model
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  # Print out the probabilities of each model being the likely explanation for the data
  #print("Probability of each model explaining the data")
  #print(strengths * 100)
  
  # Print out how many more times likely the best model is
  #print("The factor by which the best model is better than the rest")
  #print(max(strengths) / strengths)
  
  bic.data <- data.frame(cbind(scores, strengths * 100, max(strengths) / strengths, deltas))
  colnames(bic.data) <- c("scores", "probability", "factor", "deltas")
  return(bic.data)
  
  # References
  # [1] Burnham, K. P., and D. R. Anderson. 2002. Model selection and multimodel inference: a practical information-theoretic approach. 
  #     Second edition. Springer, New York, USA.
  # [2] Anderson, D. R. 2008. Model based inference in the life sciences: a primer on evidence.
  #     Springer, New York, USA.
}

# This model was tested by our primary reference paper
model.one <- function() {
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  
  # Remove any NA values from the data
  indx <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)))))
  FBN1.adipose <- FBN1.adipose[-indx]
  GLU.4wk <- GLU.4wk[-indx]
  INS.4wk <- INS.4wk[-indx]
  print(paste("Removed ", length(indx), " rows with NA values from data.", sep=""))
  
  bic.score.1 <- BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~GLU.4wk)) # Linear
  bic.score.2 <- BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~GLU.4wk+FBN1.adipose)) # Non-Li.near
  bic.score.3 <- BIC(lm(GLU.4wk~FBN1.adipose+INS.4wk)) + BIC(lm(INS.4wk~FBN1.adipose)) # Inverse Non-Linear
  bic.score.4 <- BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~FBN1.adipose)) # Independent
  
  scores <- c(bic.score.1, bic.score.2, bic.score.3, bic.score.4)
  names(scores) <- c("linear", "non.linear", "inverse.non.linear", "independent")
  
  deltas <- scores - min(scores)
  
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  bic.data <- data.frame(cbind(scores, strengths * 100, max(strengths) / strengths, deltas))
  colnames(bic.data) <- c("scores", "probability", "factor", "deltas")
  return(bic.data)
}

# This model suggests that Fbn1 affects Glucose, which affects Insulin
# Furthermore, Fbn1 also affects the expression of GPCRs
model.two <- function() {
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPR12.liver <- gene.exp("Gpr12", liver.rz)
  
  # Remove any NA values from the data
  indx.gpr12 <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPR12.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gpr12]
  GLU.4wk <- GLU.4wk[-indx.gpr12]
  INS.4wk <- INS.4wk[-indx.gpr12]
  GPR12.liver <- GPR12.liver[-indx.gpr12]
  print(paste("Removed ", length(indx.gpr12), " rows with NA values from Gpr12-related data.", sep=""))
  
  bic.score.1 <- BIC(lm(GPR12.liver~FBN1.adipose)) + BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~GLU.4wk))
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPR21.liver <- gene.exp("Rabgap1,Gpr21", liver.rz)
  
  # Remove any NA values from the data
  indx.gpr21 <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPR21.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gpr21]
  GLU.4wk <- GLU.4wk[-indx.gpr21]
  INS.4wk <- INS.4wk[-indx.gpr21]
  GPR21.liver <- GPR21.liver[-indx.gpr21]
  print(paste("Removed ", length(indx.gpr21), " rows with NA values from Gpr21-related data.", sep=""))
  
  bic.score.2 <- BIC(lm(GPR21.liver~FBN1.adipose)) + BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~GLU.4wk))  
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPRC5B.liver <- gene.exp("Gprc5b", liver.rz)
  
  # Remove any NA values from the data
  indx.gprc5b <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPRC5B.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gprc5b]
  GLU.4wk <- GLU.4wk[-indx.gprc5b]
  INS.4wk <- INS.4wk[-indx.gprc5b]
  GPRC5B.liver <- GPRC5B.liver[-indx.gprc5b]
  print(paste("Removed ", length(indx.gprc5b), " rows with NA values from Gprc5b-related data.", sep=""))
  
  bic.score.3 <- BIC(lm(GPRC5B.liver~FBN1.adipose)) + BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~GLU.4wk))
  
  scores <- c(bic.score.1, bic.score.2, bic.score.3)
  names(scores) <- c("Gpr12", "Gpr21", "Gprc5b")
  
  deltas <- scores - min(scores)
  
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  bic.data <- data.frame(cbind(scores, strengths * 100, max(strengths) / strengths, deltas))
  colnames(bic.data) <- c("scores", "probability", "factor", "deltas")
  return(bic.data)
}

# This model suggests that glucose levels influence Fbn1 expression and insulin
# Furthermore, glucose levels also influence gene expression of GPCRs
model.three <- function() {
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPR12.liver <- gene.exp("Gpr12", liver.rz)
  
  # Remove any NA values from the data
  indx.gpr12 <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPR12.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gpr12]
  GLU.4wk <- GLU.4wk[-indx.gpr12]
  INS.4wk <- INS.4wk[-indx.gpr12]
  GPR12.liver <- GPR12.liver[-indx.gpr12]
  print(paste("Removed ", length(indx.gpr12), " rows with NA values from Gpr12-related data.", sep=""))
  
  bic.score.1 <- BIC(lm(FBN1.adipose~GLU.4wk)) + BIC(lm(INS.4wk~FBN1.adipose)) + BIC(lm(GPR12.liver~FBN1.adipose))
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPR21.liver <- gene.exp("Rabgap1,Gpr21", liver.rz)
  
  # Remove any NA values from the data
  indx.gpr21 <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPR21.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gpr21]
  GLU.4wk <- GLU.4wk[-indx.gpr21]
  INS.4wk <- INS.4wk[-indx.gpr21]
  GPR21.liver <- GPR21.liver[-indx.gpr21]
  print(paste("Removed ", length(indx.gpr21), " rows with NA values from Gpr21-related data.", sep=""))
  
  bic.score.2 <- BIC(lm(FBN1.adipose~GLU.4wk)) + BIC(lm(INS.4wk~FBN1.adipose)) + BIC(lm(GPR21.liver~FBN1.adipose))
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPRC5B.liver <- gene.exp("Gprc5b", liver.rz)
  
  # Remove any NA values from the data
  indx.gprc5b <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPRC5B.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gprc5b]
  GLU.4wk <- GLU.4wk[-indx.gprc5b]
  INS.4wk <- INS.4wk[-indx.gprc5b]
  GPRC5B.liver <- GPRC5B.liver[-indx.gprc5b]
  print(paste("Removed ", length(indx.gprc5b), " rows with NA values from Gprc5b-related data.", sep=""))
  
  bic.score.3 <- BIC(lm(FBN1.adipose~GLU.4wk)) + BIC(lm(INS.4wk~FBN1.adipose)) + BIC(lm(GPRC5B.liver~FBN1.adipose))
  
  scores <- c(bic.score.1, bic.score.2, bic.score.3)
  names(scores) <- c("Gpr12", "Gpr21", "Gprc5b")
  
  deltas <- scores - min(scores)
  
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  bic.data <- data.frame(cbind(scores, strengths * 100, max(strengths) / strengths, deltas))
  colnames(bic.data) <- c("scores", "probability", "factor", "deltas")
  return(bic.data)
}

# This model compares Gprc5b between Model Two and Model Three
model.four <- function() {
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  
  # Remove any NA values from the data
  indx <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)))))
  FBN1.adipose <- FBN1.adipose[-indx]
  GLU.4wk <- GLU.4wk[-indx]
  INS.4wk <- INS.4wk[-indx]
  print(paste("Removed ", length(indx), " rows with NA values from data.", sep=""))
  
  bic.score.1 <- BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~GLU.4wk))
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPRC5B.liver <- gene.exp("Gprc5b", liver.rz)
  
  # Remove any NA values from the data
  indx.gprc5b <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPRC5B.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gprc5b]
  GLU.4wk <- GLU.4wk[-indx.gprc5b]
  INS.4wk <- INS.4wk[-indx.gprc5b]
  GPRC5B.liver <- GPRC5B.liver[-indx.gprc5b]
  
  bic.score.2 <- BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(INS.4wk~GLU.4wk)) + BIC(lm(GPRC5B.liver~FBN1.adipose))
  bic.score.3 <- BIC(lm(FBN1.adipose~GLU.4wk)) + BIC(lm(INS.4wk~FBN1.adipose)) + BIC(lm(GPRC5B.liver~GLU.4wk))
  
  scores <- c(bic.score.1, bic.score.2, bic.score.3)
  names(scores) <- c("independent", "Gprc5b.1", "Gprc5b.2")
  
  deltas <- scores - min(scores)
  
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  bic.data <- data.frame(cbind(scores, strengths * 100, max(strengths) / strengths, deltas))
  colnames(bic.data) <- c("scores", "probability", "factor", "deltas")
  return(bic.data)
}

# This model suggests that Glucose is caused by a combination of Fbn1 and GPCR expression
model.five <- function() {
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPR12.liver <- gene.exp("Gpr12", liver.rz)
  
  # Remove any NA values from the data
  indx.gpr12 <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPR12.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gpr12]
  GLU.4wk <- GLU.4wk[-indx.gpr12]
  INS.4wk <- INS.4wk[-indx.gpr12]
  GPR12.liver <- GPR12.liver[-indx.gpr12]
  
  bic.score.1 <- BIC(lm(GLU.4wk~FBN1.adipose+GPR12.liver)) + BIC(lm(INS.4wk~GLU.4wk))
  bic.score.2 <- BIC(lm(GLU.4wk~FBN1.adipose+GPR12.liver)) + BIC(lm(INS.4wk~GLU.4wk+FBN1.adipose))
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPR21.liver <- gene.exp("Rabgap1,Gpr21", liver.rz)
  
  # Remove any NA values from the data
  indx.gpr21 <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPR21.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gpr21]
  GLU.4wk <- GLU.4wk[-indx.gpr21]
  INS.4wk <- INS.4wk[-indx.gpr21]
  GPR21.liver <- GPR21.liver[-indx.gpr21]
  
  bic.score.3 <- BIC(lm(GLU.4wk~FBN1.adipose+GPR21.liver)) + BIC(lm(INS.4wk~GLU.4wk))
  bic.score.4 <- BIC(lm(GLU.4wk~FBN1.adipose+GPR21.liver)) + BIC(lm(INS.4wk~GLU.4wk+FBN1.adipose))
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPRC5B.liver <- gene.exp("Gprc5b", liver.rz)
  
  # Remove any NA values from the data
  indx.gprc5b <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPRC5B.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx.gprc5b]
  GLU.4wk <- GLU.4wk[-indx.gprc5b]
  INS.4wk <- INS.4wk[-indx.gprc5b]
  GPRC5B.liver <- GPRC5B.liver[-indx.gprc5b]
  
  bic.score.5 <- BIC(lm(GLU.4wk~FBN1.adipose+GPRC5B.liver)) + BIC(lm(INS.4wk~GLU.4wk))
  bic.score.6 <- BIC(lm(GLU.4wk~FBN1.adipose+GPRC5B.liver)) + BIC(lm(INS.4wk~GLU.4wk+FBN1.adipose))
  
  scores <- c(bic.score.1, bic.score.2, bic.score.3, bic.score.4, bic.score.5, bic.score.6)
  names(scores) <- c("Gpr12", "Gpr12.adv", "Gpr21", "Gpr21.adv", "Gprc5b", "Gprc5b.adv")
  
  deltas <- scores - min(scores)
  
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  bic.data <- data.frame(cbind(scores, strengths * 100, max(strengths) / strengths, deltas))
  colnames(bic.data) <- c("scores", "probability", "factor", "deltas")
  return(bic.data)
}

model.test <- function() {
  
  FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
  GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
  INS.4wk <- log(clinical("INS.4wk", phenotypes))
  GPRC5B.liver <- gene.exp("Gprc5b", liver.rz)
  
  # Remove any NA values from the data
  indx <- sort(unique(c(which(is.na(FBN1.adipose)), which(is.na(GLU.4wk)), which(is.na(INS.4wk)), which(is.na(GPRC5B.liver)))))
  FBN1.adipose <- FBN1.adipose[-indx]
  GLU.4wk <- GLU.4wk[-indx]
  INS.4wk <- INS.4wk[-indx]
  GPRC5B.liver <- GPRC5B.liver[-indx]
  print(paste("Removed ", length(indx), " rows with NA values from data.", sep=""))
  
  A <- FBN1.adipose
  B <- GLU.4wk
  C <- INS.4wk
  D <- GPRC5B.liver
  
  #bic.score.1 <- BIC(lm(GLU.4wk~FBN1.adipose)) + BIC(lm(GLU.4wk~GPRC5B.liver))
  #bic.score.2 <- BIC(lm(GLU.4wk~FBN1.adipose+GPRC5B.liver)) 
  #bic.score.3 <- BIC(lm(GLU.4wk~FBN1.adipose+GPRC5B.liver)) + BIC(lm(GPRC5B.liver~FBN1.adipose))
  #bic.score.4 <- BIC(lm(GLU.4wk~FBN1.adipose+GPRC5B.liver)) + BIC(lm(FBN1.adipose~GPRC5B.liver))
  
  bic.score.1 <- BIC(lm(B~A)) + BIC(lm(B~D)) + BIC(lm(C~B))
  #bic.score.2 <- BIC(lm(B~A+D))
  bic.score.2 <- BIC(lm(B~A+D)) + BIC(lm(C~B))
  
  scores <- c(bic.score.1, bic.score.2)
  names(scores) <- c("independent", "linear")
  
  deltas <- scores - min(scores)
  
  strengths <- exp(-0.5 * deltas) / sum(exp(-0.5 * deltas))
  
  bic.data <- data.frame(cbind(scores, strengths * 100, max(strengths) / strengths, deltas))
  colnames(bic.data) <- c("scores", "probability", "factor", "deltas")
  return(bic.data)
}

round(model.one(), digits=2)
round(model.two(), digits=2)
round(model.three(), digits=2)
round(model.four(), digits=2)
round(model.five(), digits=2)