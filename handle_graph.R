rm(list=ls())

require(igraph)

pplot <- function(g,cc=0,...) plot(g,vertex.color=cc,...)
p <- 20

set.seed(123)
g <- erdos.renyi.game(p,3/p)
V(g)$name <- letters[1:p]
vec <- sample(1:p,0.5*p)
col <- ( (1:p) %in% vec ) * 2

a <- induced.subgraph(g,v=vec)
b <- induced.subgraph(g,v=(1:p)[-vec])


layout(matrix(c(1,1,2,3),2,2))
pplot(g,cc=col)
pplot(a)
pplot(b,cc=2)

A <- get.adjacency(a)
indchar <- apply(which(A==0 & upper.tri(A),arr.ind=T),2,function(x) V(a)$name[x])
foo <- function(ll) lapply(ll,function(o) V(g)$name[o[-c(1,length(o))]])
ulist1 <- apply(indchar,1,function(v) foo(get.all.shortest.paths(g,v[1],v[2])$res ))
ulist2 <- apply(indchar,1,function(v) (get.all.shortest.paths(g,v[1],v[2])$res ))

bool <- sapply(ulist1,function(ll) all(sapply(ll,function(v) all(v %in% V(b)$name) )))

nbe <- length(E(a)) ; coledge <- rep(0,nbe)
a[from=indchar[bool,1],to=indchar[bool,2]] <- 1
coledge <- c(coledge,rep(1,length(E(a))-nbe))

x11()
pplot(a,edge.color=coledge+1)

x11()
pplot(g,cc=col)


