###########################################################
#  This is a script to demo the use of Attie BTBR eQTL data
#  This script looks at lipid phenotypes
#  updated January 7, 2014 - GAC
###########################################################
# 
source("src/configuration.R")

names(phenotypes.rz)
#we are going to look at LDL, HDL and CHOL

# IMPORTANT!! SET UP THE DATA HERE
#set up the cross object with out phenotypes
summary(f2g)
names(f2g$pheno)
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],phenotypes.rz[,c("LDL","HDL","CHOL")])
names(f2g$pheno)

#look at histrograms of data by sex
graphics.open()
qplot(CHOL, facets=Sex~., data=f2g$pheno)
graphics.open()
qplot(HDL, facets=Sex~., data=f2g$pheno)
graphics.open()
qplot(LDL, facets=Sex~., data=f2g$pheno)
# note higher CHOL and HDL in males, also LDL but less so

#look at raw data cholesterol traits
#log scaling is helpful
graphics.open()
pairs(log(phenotypes[,c("LDL","HDL","CHOL")]),upper.panel=panel.cor,diag.panel=panel.hist)
#normally LDL is computed trait
graphics.open()
qplot(CHOL-HDL, LDL, data=phenotypes)
graphics.open()
qplot(LDL+HDL, CHOL, data=phenotypes)		
#LDL appears to be directly measured here
#ask Mark about this

#look at transformed cholesterol traits
graphics.open()
pairs(f2g$pheno[,c("LDL","HDL","CHOL")],	upper.panel=panel.cor,diag.panel=panel.hist)
#all three are positively correlated
#HDL and CHOL are tightly correlated 

##########################################################
#compute principle components	
pc.chol <- princomp(~CHOL+HDL+LDL, f2g$pheno, na.action=na.exclude,cor=TRUE)
names(pc.chol)
summary(pc.chol)

#look at loadings to interpret the components
pc.chol$loadings

#Loadings:
#     Comp.1 Comp.2 Comp.3
#CHOL -0.595  0.391  0.702
#HDL  -0.597  0.371 -0.712
#LDL  -0.539 -0.842 

# comp 1 is a weighted average about equal weights, note negative sign
# comp 2 contrasts HDL&CHOL to LDL, this is our "cholesterol ratio"
# comp 3 is the difference HDL vs CHOL - a computed LDL

#default plot is "scree" tells us how important each component is
graphics.open()
plot(pc.chol)

#biplots show us the data in rotated coordinates
graphics.open()
biplot(pc.chol)
graphics.open()
biplot(pc.chol,choices=c(1,3))
graphics.open()
biplot(pc.chol,choices=c(2,3))

#add PCs to our pheno data
f2g$pheno = transform(f2g$pheno, PC1 = pc.chol$scores[,1], 
	PC2 = pc.chol$scores[,2], PC3 = pc.chol$scores[,3])

##########################################################
# genome scans
# set up
f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

#convenient to keep "sex" handy as a numeric variable
sex <- as.numeric(f2g$pheno$Sex)-1

#run permutations -- may take a long time.  Ignore warnings
# running only for one trait - lazy
f2g.perm1a <-scanone(f2g, pheno.col="CHOL", addcovar=sex,method="hk",n.perm=100,perm.Xsp=TRUE)
# save results - so you can reload the permutation results later 
save(list="f2g.perm1a", file="f2g_perm1a.Rdata")
#load(file="f2g_perm1a.Rdata")

#scan all traits with sex as an additive covariate
f2g.scan1a <- scanone(f2g, pheno.col=4:6, addcovar=sex, method="hk")

#plot the genome scans
graphics.open()
par(mfrow=c(3,1))
plot(f2g.scan1a,lodcolumn=1)
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1a,lodcolumn=2)
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1a,lodcolumn=3)
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
###
graphics.open()
par(mfrow=c(3,1))
plot(f2g.scan1a,lodcolumn=4)
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1a,lodcolumn=5)
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1a,lodcolumn=6)
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1a, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")


#tabular summary	
summary(f2g.scan1a)
summary(f2g.scan1a, perms=f2g.perm1a, alpha=0.10)
summary(f2g.scan1a, perms=f2g.perm1a, alpha=0.10, format="tabByCol",ci.function="lodint")
#
#confidence interval on Chr 17
lodint(f2g.scan1a,17,lodcolumn=3)
#          
#rs8256996  17 11.82007 
#rs3700924  17 13.41000 
#rs6358703  17 14.61075 

#effect plot
graphics.open()
plotPXG(f2g, find.marker(f2g,17,13.41), pheno.col="CHOL")
# BTBR is recessive low  == B6 is dominant high

#look at other peaks for CHOL and break out by sex
graphics.open()
par(mfrow=c(1,3))
effectplot(f2g, pheno.col="CHOL", mname1="Sex",mark1=sex,
	mname2 = find.marker(f2g,7,19.14),ylim=c(-0.9,0.7))
effectplot(f2g, pheno.col="CHOL", mname1 ="Sex",mark1=sex,
	mname2= find.marker(f2g,12,4.86),ylim=c(-0.9,0.7))
effectplot(f2g, pheno.col="CHOL", mname1 ="Sex",mark1=sex,
	mname2= find.marker(f2g,17,13.41),ylim=c(-0.9,0.7))
	
##########################################################
# scan with sex as an interactive covariate

#run permutations
#f2g.perm1b <- scanone(f2g, pheno.col="CHOL", intcovar=sex, method="hk", 
#	n.perm=100,perm.Xsp=TRUE)	
# save results 
#save(list="f2g.perm1b", file="f2g_perm1b.Rdata")
load(file="f2g_perm1b.Rdata")

# run scans
f2g.scan1b <- scanone(f2g, pheno.col=4:9, intcovar=sex, method="hk")
#

#plot the genome scans
graphics.open()
par(mfrow=c(3,1))
plot(f2g.scan1b, lodcolumn=1)
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1b, lodcolumn=2)
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1b, lodcolumn=3)
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
graphics.open()
par(mfrow=c(3,1))
plot(f2g.scan1b, lodcolumn=4)
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1b, lodcolumn=5)
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#
plot(f2g.scan1b, lodcolumn=6)
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1b, perms=f2g.perm1a, alpha=0.63,lty="dashed",lwd=2,col="green")
#

#
summary(f2g.scan1b)
summary(f2g.scan1b, perms=f2g.perm1b, alpha=0.10)
summary(f2g.scan1b, perms=f2g.perm1b, alpha=0.10, format="tabByCol",ci.function="lodint")
#

#effect plot
quartz(width=10, height=5)
par(mfrow=c(1,2))
effectplot(f2g, pheno.col="CHOL", mname1="Sex", mname2=find.marker(f2g,8,58.28),
	main = "Interaction of Sex and Chr 8")
#
effectplot(f2g, pheno.col="CHOL", mname1="Sex", mname2= find.marker(f2g,17,14.31),
	main = "Interaction of Sex and Chr 17")

##########################################################
#pairscan with sex as an additive covariate

# setup for pair scan at a lower density = 2cM step
f2g <- calc.genoprob(f2g, step=2, stepwidth="fixed", map.function="c-f", err=0.002)

#pairscan with sex as an additive covariate
#f2g.scan2a <- scantwo(f2g, pheno.col="CHOL", addcovar=sex, method="hk")

#save(list="f2g.scan2a", file="f2g_scan2a.Rdata")
load(file="f2g_scan2a.Rdata")

#plot
graphics.open()
plot(f2g.scan2a)
graphics.open()
plot(f2g.scan2a,chr=c(5,7,12,17))
graphics.open()
plot(f2g.scan2a,chr=c(5,7,12,17),upper="cond-int",lower="cond-add")

#report
summary(f2g.scan2a, what="best", thresholds=c(9.1, 7.1, 6.3, 6.3, 3.3))

summary(f2g.scan2a, allpairs=FALSE)

#set up for effect plots
f2g <- sim.geno(f2g,step=2)

# effect plot
graphics.open()
par(mfrow=c(2,2))
effectplot(f2g, pheno.col="CHOL", mname1=find.marker(f2g,5,34.26),
	ylim=c(-1.4,0.8),main="",xlab="Chr 5 @ 34cM")
#
effectplot(f2g, pheno.col="CHOL", mname2=find.marker(f2g,5,34.26), 
	mname1= find.marker(f2g,17,13.3),ylim=c(-1.4,0.8),
	main="",xlab="Chr 5 @ 34cM",legend.lab="Chr 17",add.legend=FALSE)
#
effectplot(f2g, pheno.col="CHOL", mname1=find.marker(f2g,17,13.3),
	ylim=c(-1.4,0.8),main="",xlab="Chr 17 @ 13cM")
#
effectplot(f2g, pheno.col="CHOL", mname1=find.marker(f2g,5,34.26), 
	mname2=find.marker(f2g,17,13.3),ylim=c(-1.4,0.8),
	main="",xlab="Chr 17 @ 13cM",legend.lab="Chr 5",add.legend=TRUE)
	

##########################################################
# conditional scan with chrom 17 QTL as a covariate

#find the genoprobs
Q17 <- f2g$geno[[17]]$prob[,find.marker(f2g,17,13.3),]

#create additive covariate
ac <- cbind(sex, Q17[,3]) 

#run permutations -- may take a long time.  Ignore warnings
#f2g.perm1c <-scanone(f2g, pheno.col="CHOL", addcovar=ac, method="hk",
#	n.perm=100,perm.Xsp=TRUE)
# save results - so you can reload the permutation results later 
#save(list="f2g.perm1c", file="f2g_perm1c.Rdata")
load(file="f2g_perm1c.Rdata")

# conditional scan with sex and Chr 17 as an additive covariate
f2g.scan1c <- scanone(f2g, pheno.col="CHOL", addcovar=ac, method="hk")

#plot the genome scan
graphics.open()
plot(f2g.scan1c)
add.threshold(f2g.scan1c, perms=f2g.perm1c, 
	alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(f2g.scan1c, perms=f2g.perm1c, 
	alpha=0.63,lty="dashed",lwd=2,col="green")

#tabular summary 	
summary(f2g.scan1c, perms=f2g.perm1c, alpha=0.10)

#multiple QTL models
qtls <- makeqtl(f2g, chr=c(5,7,12,17), 
	pos=c(35.45, 20.14, 4.77, 13.3))
	
fit1 <- fitqtl(f2g,pheno.col="CHOL",qtls,covar=as.data.frame(sex),
	y~sex+Q1+Q2+Q3+Q4)
summary(fit1)

###############################################################
# 
# what genes are in the chr 17 CI?   

pmap[[17]]["rs8256996"]
pmap[[17]]["rs6358703"]

tmp <- subset(annot, select=c("a_gene_id","gene_symbol","chromosome","start","end"))
tmp2 <- subset(tmp, chromosome=="chr17"  &
	end>pmap[[17]]["rs8256996"]*1000000 & start<pmap[[17]]["rs6358703"]*1000000)







