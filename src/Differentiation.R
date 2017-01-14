#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: November 20 2016                      ##
#################################################

source("src/configuration.R")
source("src/BICScores.R")

names(f2g$pheno)
f2g$pheno <- f2g$pheno[,c("MouseNum", "Sex", "pgm")]
ids <- seq(nrow(f2g$pheno))
FBN1.adipose <- gene.exp("Fbn1", adipose.rz)
GLU.4wk <- log(clinical("GLU.4wk", phenotypes))
f2g$pheno <- cbind(ids, f2g$pheno, FBN1.adipose, GLU.4wk)
names(f2g$pheno)

ids.selected <- f2g$pheno$ids[which(f2g$pheno$GLU.4wk>mean(GLU.4wk))]
f2g$pheno <- f2g$pheno[ids.selected,]

for (i in seq(19)) {
  f2g$geno[[i]]$data <- f2g$geno[[i]]$data[ids.selected,]
}
f2g$geno$X$data <- f2g$geno$X$data[ids.selected,]

f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

sex <- as.numeric(f2g$pheno$Sex) - 1

FBN1.adipose.scan1 <- scanone(f2g, pheno.col=5:length(f2g$pheno), addcovar=sex, method="hk")
FBN1.adipose.perms <- scanone(f2g, pheno.col=5, addcovar=sex, method="hk", n.perm=100, perm.Xsp=TRUE)

summary(FBN1.adipose.scan1, format="tabByChr", perms=FBN1.adipose.perms, alpha=0.3)

# Chr9@47.3