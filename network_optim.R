rm(list=ls())

.libPaths("/volatile/R/x86_64-pc-linux-gnu-library/2.15/")
setwd("/volatile/OverOptim")

load("EColi/EColi.RData")
#ind <- sample(nrow(x),100)
#X <- t(as.matrix(x)[ind,])
#true.pcor <- adjmat[ind,ind]
X <- t(as.matrix(x))
true.pcor <- adjmat

 ###########
## GeneNet ##
 ###########

# load GeneNet library
library("GeneNet")

pc <- ggm.estimate.pcor(X)
x.edges <- network.test.edges(pc, direct=TRUE, fdr=TRUE)
x.net <- extract.network(x.edges)
dim(x.net)
res.genenet <- matrix(0,ncol(X),ncol(X))
res.genenet[as.matrix(x.net[,c("node1","node2")])] <- 1

 ##########
## PC Alg ##
 ##########

require(pcalg)
## define independence test (partial correlations)
indepTest <- gaussCItest
## define sufficient statistics
suffStat <- list(C = cor(X), n = nrow(X))
## estimate CPDAG
alpha <- 0.05
pc.fit <- pc(suffStat, indepTest, p=ncol(X), alpha, verbose = TRUE)
res.pcalg <- as(pc.fit@graph,"matrix")

 ##############
## Comparison ##
 ##############

require(parcor)
uptrind <- upper.tri(true.pcor) 
pc0 <- (true.pcor[uptrind]!=0) + 0
pc1 <- (res.genenet[uptrind]!=0) + 0
pc2 <- (res.pcalg[uptrind]!=0) + 0

n00 <- sum( (pc0+pc1+pc2)==0)
n11 <- sum( (pc0+pc1+pc2)==3)
n01 <- sum( (pc0-pc1)==1 & (pc0+pc2)==2 )
n10 <- sum( (pc0+pc1)==2 & (pc0-pc2)==1 )

(n00*n11)/(n10*n01)
(n00/n10)*(n11/n01)
final.tab <- matrix(as.numeric(c(n00,n10,n01,n11)),2,2)

ntot <- sum(final.tab)/100

big.final.tab <- rbind(cbind(final.tab/ntot,rowSums(final.tab)/ntot),c(colSums(final.tab)/ntot,100))

mat <- round(final.tab/ntot*100)


compute.prob <- function(mat) {
 N <- sum(mat)
 print(dim(mat))
 n10 <- mat[1,1] + mat[1,2]
 n01 <- mat[1,1] + mat[2,1]
 n02 <- mat[1,2] + mat[2,2]
 n20 <- mat[2,1] + mat[2,2]
 theta <- mat[1,1]*mat[2,2]/(mat[1,2]*mat[2,1])
 X <- 1:N
 num <- choose(n10,n01-X)*choose(n20,X)*theta^(n01-X) 
 den <- sum( sapply(max(0,n01-n20):min(n01,n10), function(u) 
    choose(n10,u)*choose(n20,n01-u)*theta^(u)) )
 return(sum(num/den))
}

compute.prob(round(final.tab/ntot))






 ######################################
## Apply GeneNet to the whole dataset ##
##      and a subsample of it         ##
 ######################################

#ggm.wrap <- function(X){
 # pc <- ggm.estimate.pcor(X,verbose=F)
 # x.edges <- network.test.edges(pc, direct=TRUE, fdr=TRUE,verbose=F,plot=F)
 # x.net <- extract.network(x.edges, verbose=F)
 # res.genenet <- matrix(0,ncol(X),ncol(X))
 # res.genenet[as.matrix(x.net[,c("node1","node2")])] <- 1
 # return(res.genenet)
#}


ggm.wrap <- function(X){
  pc <- ggm.estimate.pcor(X,verbose=F)
  cutoff <- fdrtool(pc[upper.tri(pc)], statistic = "correlation",
		    verbose=F,plot=F)$param[1]
  return(adja = (abs(pc) > cutoff) + 0)
}



# load GeneNet library
require(GeneNet)

B <- 100
N <- nrow(X)
pctot <- pcsub <- NULL

res.cvnet <- NULL

pctot <- ggm.wrap(X)
uptrind <- upper.tri(true.pcor) 
v0 <- (true.pcor[uptrind] != 0) +0
v1 <- (pctot[uptrind] != 0) +0

pb <- txtProgressBar(style=3)
for (i in 1:B) {
  ind <- sample(N,0.9*N)
  pcsub <- ggm.wrap(X[ind,])
  
  v2 <- (pcsub[uptrind] != 0) +0
  n00 <- sum( (v0+v1+v2)==0)
  n11 <- sum( (v0+v1+v2)==3)
  n01 <- sum( (v0-v1)==1 & (v0+v2)==2 )
  n10 <- sum( (v0+v1)==2 & (v0-v2)==1 )
  
  res.cvnet[[i]] <- matrix(c(n00,n10,n01,n11),2,2)
  
  setTxtProgressBar(pb, i/B)
}
close(pb)

computetheta <- function(M) (M[1,1]*M[2,2])/(M[1,2]*M[2,1])

final.tab <- Reduce("+",res.cvnet)/100
thetahat <- computetheta(final.tab)

ntot <- sum(final.tab)/100

big.final.tab <- rbind(cbind(final.tab/ntot,rowSums(final.tab)/ntot),c(colSums(final.tab)/ntot,100))


prob.hat <- compute.prob(final.tab/ntot)




























require(parcor)






require(pacose)
require(GeneNet)
# 
ind <- sample(nrow(x),100)
X <- t(as.matrix(x)[ind,])
true.pcor <- adjmat[ind,ind]


# pls.net
pls.res <- pls.net(X,verbose=TRUE,k=3)
ridge.res <- ridge.net(X,verbose=TRUE,k=3,cv.method="HKB")
#ridge.res <- ridge.net(X,verbose=TRUE,k=2,cv.method="CV")
# ridge.net

# lasso.net

# aracne

perf1 <- performance.pcor(pls.res$pcor, true.pcor, fdr=TRUE, verbose=FALSE, plot=FALSE)
perf2 <- performance.pcor(ridge.res$pcor, true.pcor, fdr=TRUE, verbose=FALSE, plot=FALSE)

uptrind <- upper.tri(true.pcor) 
pc0 <- (true.pcor[uptrind]==0) +0
pc1 <- (ridge.res$pcor[uptrind]==0) +0
pc2 <- (pls.res$pcor[uptrind]==0) +0


n00 <- sum( (pc0+pc1+pc2)==0)
n11 <- sum( (pc0+pc1+pc2)==3)
n01 <- sum( (pc0-pc1)==1 & (pc0+pc2)==2 )
n10 <- sum( (pc0+pc1)==2 & (pc0-pc2)==1 )

(n00*n11)/(n10*n01)


