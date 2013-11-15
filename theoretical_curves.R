###############################################
### Plot the probability distribution of the ##
### number of new interactions found with    ##
### several different sets of parameters.    ##
###############################################

rm(list=ls())

N <- 50
pp <- 0.5 # the percentage of TRUE discoreved edges

# Compute the probabilities
proba <- function(x,theta) {
  n10 <- pp*N
  n01 <- pp*N
  n20 <- (1-pp)*N
  n02 <- (1-pp)*N
  num <- choose(n10,n01-x)*choose(n20,x)*theta^(n01-x) 
  den <- sum( sapply(max(0,n01-n20):min(n01,n10), function(u) 
    choose(n10,u)*choose(n20,n01-u)*theta^(u)) )
  return(num/den)
}

x <- 0:N
thetavec <- seq(0.01,100,l=1000)#seq(0.01,10,l=1000)


# Not presented in the article: several distributions
pdf("DistribTheta.pdf")
par(family="serif")
sapply(thetavec,function(theta) plot(x,proba(x,theta),type="h",#ylim=c(0,.25),
                                            xlab=expression(italic(x)),
                                            ylab=expression(P(n[21]==x)),
                                            main=bquote(theta==.(theta) )))
dev.off()
     
# P(n12>0) as a function of theta
pdf("PxSup0.pdf",width=15,height=10)
par(family="serif",cex=2.2)
probs <- sapply(thetavec,function(theta) sum(proba(1:N,theta)) )
plot(thetavec,(probs),type="l",xlab=expression(theta),#log="y",
     ylab=expression(P(n[21]>0)) ,lwd=2)
dev.off()

xproba <- function(x,theta) x*proba(x,theta)

# Mean as a function of theta
pdf("mean.pdf",width=15,height=10)
moy <- sapply(thetavec,function(theta,x) sum(xproba(x,theta)),0:N)
par(family="serif",cex=2.2)
plot(thetavec,moy,type="l",xlab=expression(theta),ylab="Mean",lwd=2)
dev.off()






# ######### TO CLEAN ###########
# getmedi <- function(i,mm,tt) {
#   Fx <- mm[,i]
#   xprobainf <- xproba( max((0:N)[Fx <=0.5]), tt[i])
#   xprobasup <- xproba( min((0:N)[Fx > 0.5]), tt[i])
#   return(xprobainf+xprobasup)
# }
# #mat <- sapply(thetavec,function(theta) cumsum(proba(0:N,theta)) )
# #medi <- sapply(1:ncol(mat), getmedi ,mat, thetavec)
# 
# pdf("manymeans.pdf")
# NN <- floor(seq(30,100,l=10))
# thetavec <- seq(0.01,5,l=500)
# moy <- sapply(thetavec,function(theta,x) sum(xproba(x,theta)),0:N)
# par(family="serif")
# plot(thetavec,moy,type="n",ylim=c(0,30),xlab=expression(theta),ylab="Mean")
# for(N in NN) {
#   pp <- 30/N
#   
#   proba <- function(x,theta) {
#     n10 <- pp*N
#     n01 <- pp*N
#     n20 <- (1-pp)*N
#     n02 <- (1-pp)*N
#     num <- choose(n10,n01-x)*choose(n20,x)*theta^(n01-x) 
#     den <- sum( sapply(max(0,n01-n20):min(n01,n10), function(u) 
#       choose(n10,u)*choose(n20,n01-u)*theta^(u)) )
#     return(num/den)
#   }
#   x <- 0:N
#   moy <- sapply(thetavec,function(theta,x) sum(xproba(x,theta)),0:N)
#   lines(thetavec,moy,type="l",col=N)
# }
# dev.off()
