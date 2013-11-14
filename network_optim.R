load("../data/EColi/EColi.RData")

require(parcor)

X <- t(as.matrix(x))

# pls.net
pc<-pls.net(X,ncomp=10,k=5)

# ridge.net

# lasso.net

# aracne

