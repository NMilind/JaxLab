
########################################################
# functions to jazz up the pairs() all pairwise scatterplots
# some useful plotting functions
# see documentation for "pairs" function 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use="complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor*0.5, col=c("gray60", "black")[(abs(r)>0.65)+1])
}
#
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2],0,1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


#########################################################
# correlation "scan" functions
# x is a matrix or dataframe with lots varaibles 
# y is a vector 
# cov is a covariate 

# compute correlation of y with each element of x
mycorr <- function(x,y){
	apply(x,2,"cor",y,use="complete.obs")
}

# compute residuals to adjust x wrt covariate
myresids <- function(x, cov){
	residuals(lm(x~cov,na.action=na.exclude))
}

# compute correlation of y with each element of x after adjusting for covariate cov
mycorr.adj <- function(x,y,cov){
	x.adj <- apply(x,2,"myresids",cov)
	y.adj <- myresids(y,cov)
	apply(x.adj,2,"cor",y.adj,use="complete.obs")
}

# a FWER adjustment permutation test for correlation scans
mycorr.permute <- function(x,y,n){
	max.cor <- NULL
	for(i in 1:n){
			max.cor <- c(max.cor, max(abs(mycorr(x,sample(y))),na.rm=TRUE))
	}
	max.cor
}

# a FWER adjustment permutation test for adjusted correlation scans
mycorr.adj.permute <- function(x,y,cov,n){
	x.adj <- apply(x,2,"myresids",cov)
	y.adj <- myresids(y,cov)
	max.cor <- NULL
	for(i in 1:n){
			max.cor <- c(max.cor, max(abs(mycorr(x.adj,sample(y.adj))),na.rm=TRUE))
	}
	max.cor
}

R2Z <- function(r) log((1+r)/(1-r))

# plot a histogram of data with normal distribution
histnorm <- function(x) {
        lab <- deparse(substitute(x))
        h<-hist(x, breaks=10, col="red", main=paste("Histogram of ", lab))
        xfit<-seq(min(x),max(x),length=40) 
        yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
        yfit <- yfit*diff(h$mids[1:2])*length(x) 
        lines(xfit, yfit, col="blue", lwd=2)
}


### Define functions that we will use in the analysis below

####
# some useful plotting functions
# see documentation for "pairs" function
# note that you can adjust the multiplier on cex.cor
#    and the threshold for gray vs black
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*1.0, col=c("gray60", "black")[(abs(r)>0.5)+1])
}
#
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

####
# fit causal models to a triplet with BIC scoring
# X is a transcript used here as first argument to make "apply" easy
# Y is a clincal
# Q is a genotype (factor or numeric)

triple.fit <- function(X,Y,Q){
  #remove any rows with missing values
  indx <- sort(unique(
    c(which(is.na(X)),which(is.na(Y)),which(is.na(Q)))
  ))
  X <- X[-indx]
  Y <- Y[-indx]
  Q <- Q[-indx]
  
  # fit models and compute scores
  b1 <- BIC(lm(X~Q)) + BIC(lm(Y~Q))	#independent X<-Q->Y
  b2 <- BIC(lm(X~Y)) + BIC(lm(Y~Q))	#reactive	 Q->Y->X
  b3 <- BIC(lm(X~Q)) + BIC(lm(Y~X))	#causal		 Q->X->Y
  b4 <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X)) #complex
  scores <- c(b1,b2,b3,b4)
  names(scores) <- c("independent","reactive","causal","complex")
  scores
}

gene.data.exists <- function(symbol) {
  
  return(length(annot$gene_symbol[which(annot$gene_symbol==symbol)]) > 0)
}