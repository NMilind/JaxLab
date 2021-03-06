#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

runQTLNet <- function(genes.adipose=c(), genes.liver=c(), genes.hypo=c(), genes.islet=c(), phenos=c(), gene.chrs=c(), gene.pos=c()) {
  
  #################################################
  ## ENVIRONMENTAL SETUP                         ##
  #################################################
  
  # Clear memory for QTLNet usage
  gc()
  
  #################################################
  ## SETUP GENES AND PHENOTYPES                  ##
  #################################################
  
  f2g$pheno <- f2g$pheno[,c("MouseNum","Sex","pgm")]
  for (i in 1:length(genes.adipose)) {
    f2g$pheno <- cbind(f2g$pheno, adipose.rz[,annot$a_gene_id[which(annot$gene_symbol==genes.adipose[i])]])
  }
  for (i in 1:length(genes.liver)) {
    f2g$pheno <- cbind(f2g$pheno, liver.rz[,annot$a_gene_id[which(annot$gene_symbol==genes.liver[i])]])
  }
  for (i in 1:length(genes.hypo)) {
    f2g$pheno <- cbind(f2g$pheno, hypo.rz[,annot$a_gene_id[which(annot$gene_symbol==genes.hypo[i])]])
  }
  for (i in 1:length(genes.islet)) {
    f2g$pheno <- cbind(f2g$pheno, islet.rz[,annot$a_gene_id[which(annot$gene_symbol==genes.islet[i])]])
  }
  f2g$pheno <- cbind(f2g$pheno, phenotypes.rz[,phenos])
  if (length(gene.chrs) != 0) {
    for (i in 1:length(gene.chrs)) {
      f2g$pheno <- cbind(f2g$pheno, genotype(chr=gene.chrs[i], pos=gene.pos[i]))
    }
  }
  vecs = paste(paste(rep("Chr", length(gene.chrs)), gene.chrs, sep=""), gene.pos, sep="Pos")
  names(f2g$pheno) <- c(c("MouseNum", "Sex", "pgm"), genes.adipose, genes.liver, genes.hypo, genes.islet, phenos, vecs)
  names(f2g$pheno)
  
  #################################################
  ## SETUP WORKING DIRECTORY                     ##
  #################################################
  
  do.call(file.remove, list(list.files("data/net-data/", full.names=TRUE)))
  
  #################################################
  ## SETUP FOR NETWORK                           ##
  #################################################
  
  # Clear memory again
  gc()
  # Range of phenotypes that are to be used
  pheno.col <- 4:ncol(f2g$pheno)
  # Maximum number of parents
  max.parents <- ncol(f2g$pheno) - 4
  # Number of possible scanone operations that will be performed
  size.qtlnet(pheno.col, max.parents)
  # Creates a list of parent sets that will be covariates
  parents <- parents.qtlnet(pheno.col, max.parents)
  # Groups the parents into approximately equal sizes for parallel computation
  groups <- group.qtlnet(parents=parents, group.size=8)
  # Save the data
  save(f2g, pheno.col, max.parents, parents, groups, file="data/net-data/Step1.RData", compress=TRUE)
  pa <- summary(parents)
  
  #################################################
  ## RUN SCANONE ON EACH PHENOTYPE               ##
  #################################################
  
  # Precompute BIC Scores using scanone and parallel computing
  # Create empty array
  N <- rep(NA, nrow(groups))
  for (i in 1:nrow(groups))
    N[i] <- sum(pa[seq(groups[i,1], groups[i,2]),2])
  
  # Calculate the BIC score for each row in groups
  for (i in seq(nrow(groups))) {
    bic <- bic.qtlnet(f2g, 
                      pheno.col, 
                      threshold=3.04, 
                      max.parents=max.parents, 
                      parents=parents[seq(groups[i,1], groups[i,2])], 
                      verbose=FALSE)
    cat("BIC Groups Row", i, "\n")
    save(bic, file=paste("data/net-data/bic", i, ".RData", sep=""), compress=TRUE)
  }
  
  #################################################
  ## PARALLELIZED MARKOV CHAIN MODEL             ##
  #################################################
  
  load(file="data/net-data/Step1.RData")
  
  # Read and save the BIC scores into one object
  bic.group <- list()
  for (i in seq(nrow(groups))) {
    load(file=paste("data/net-data/bic", i, ".RData", sep=""))
    bic.group[[i]] <- bic
    cat("Compiled Group", i, "\n")
  }
  
  saved.scores <- bic.join(f2g, pheno.col = 4:ncol(f2g$pheno), bic.group, max.parents = ncol(f2g$pheno)-4)
  print(saved.scores)
  
  # Generate the Markov Chain Model
  set.seed(54321)
  # Number of parallel runs
  n.runs <- 5
  
  for (i in seq(n.runs)) {
    cat("Markov Chain Run", i, "\n")
    mcmc <- mcmc.qtlnet(f2g, 
                        pheno.col=4:ncol(f2g$pheno), 
                        threshold=3.04, 
                        thinning=1, 
                        max.parents=max.parents, 
                        saved.scores=saved.scores, 
                        init.edges=NULL)
    save(mcmc, file=paste("data/net-data/mcmc", i, ".RData", sep=""), compress=TRUE)
  }
  
  #################################################
  ## POST-PROCESSING                             ##
  #################################################
  
  outs.qtlnet <-list()
  for (i in seq(n.runs)){
    load(paste("data/net-data/mcmc", i, ".RData", sep=""))
    outs.qtlnet[[i]] <- mcmc
  }
  
  # Concatenate Outputs from Various Runs
  output <- c.qtlnet(outs.qtlnet) 
  
  # Summary Table of Causal Directionality
  print(summary(output))
  # Plot the Causal Network
  plot(output)
  
  # Save for Future Use
  save(file="data/net-data/output.RData", output)
}

#################################################
## LOAD OUTPUT AND VIEW                        ##
#################################################

loadQTLOutput <- function(loadfile="data/net-data/output.RData") {
  
  load(file=loadfile)
  print(summary(output))
  plot(output)
}

#runQTLNet(genes.liver=c("Gprc5b"), phenos=c("HOMA.10wk"), gene.chrs=c(2,3), gene.pos=c(61.38,49.5))
#runQTLNet(genes.liver=c("Gpr21"), phenos=c("HOMA.10wk"), gene.chrs=c(2,3), gene.pos=c(61.38,49.5))
#runQTLNet(genes.adipose=c("Amy1", "Fbn1"), phenos=c("INS.10wk"))
