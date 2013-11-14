rm(list=ls())

# setwd("/volatile/PACOSE/")
# setwd("F:/Network_overoptimization_VG_VJ_ALB/OverOptim/")

require(gdata)
tmptab <- read.table("E-GEOD-29536.sdrf.csv",header=TRUE,sep="\t",as.is=TRUE)[,c("Source.Name","FactorValue..STATUS.")]
tmptab[,1] <- gsub(" 1","",tmptab[,1])
names(tmptab) <- c("ID","Status")
tab <- subset(tmptab,grepl("[Hh]ealthy|[Rr]ecovery",Status))
filenames <- sapply(tab$ID,grep,dir("GEO-29536",full.names=TRUE),value=TRUE)

xlist <- lapply(filenames,function(ff) read.table(ff,header=TRUE))

# check the probes are the same for each array
probeslist <- sapply(xlist,rownames)
res <- TRUE
for (i in 2:116) res+(all(probeslist[,1]==probeslist[,i]))
res

# build matrix of expression values
Xraw <- t(sapply(xlist,function(mat) mat[,1]))

# build matrix of significativity for the expression values
pvals <- t(sapply(xlist,function(mat) mat[,2]))

colnames(Xraw) <- colnames(pvals) <- rownames(xlist[[1]])
rownames(Xraw) <- rownames(pvals) <- tab$ID

X <- scale(log(Xraw[,apply(pvals,2,function(vec) all(vec < 0.05))]))

plot(prcomp(X,scale=TRUE)$x[,1:2],col=as.numeric(factor(tab$Status)),pch=16)

jpeg("tmp.jpg")
boxplot(t(X),col=as.numeric(factor(tab$Status)))
dev.off()

## Annotations + network extraction

require(illuminaHumanv2.db)
require(KEGGgraph)

pbid <- unlist(as.list(illuminaHumanv2SYMBOL)[colnames(X)])
kegg2pid <- unlist(as.list(illuminaHumanv2PATH2PROBE))

f1 <- "/media/E832-5FD7/xml_2/hsa04612.xml"
f2 <- "/media/E832-5FD7/xml_2/hsa04660.xml"
f3 <- "/media/E832-5FD7/xml_2/hsa04662.xml"
f4 <- "/media/E832-5FD7/xml_2/hsa04666.xml"
f5 <- "/media/E832-5FD7/xml_2/hsa04670.xml"
f6 <- "/media/E832-5FD7/xml_2/hsa04810.xml"

f1 <- "xml_2/hsa04612.xml"
f2 <- "xml_2/hsa04660.xml"
f3 <- "xml_2/hsa04662.xml"
f4 <- "xml_2/hsa04666.xml"
f5 <- "xml_2/hsa04670.xml"
f6 <- "xml_2/hsa04810.xml"

g1 <- parseKGML2Graph(f1,genesOnly=TRUE)
g2 <- parseKGML2Graph(f2,genesOnly=TRUE)
g3 <- parseKGML2Graph(f3,genesOnly=TRUE)
g4 <- parseKGML2Graph(f4,genesOnly=TRUE)
g5 <- parseKGML2Graph(f5,genesOnly=TRUE)
g6 <- parseKGML2Graph(f6,genesOnly=TRUE)
glist <- list(hsa04612=g1,hsa04660=g2,hsa04662=g3,hsa04666=g4,hsa04670=g5,hsa04810=g6)


V <- unique(as.vector(unlist(sapply(glist, nodes))))
E <- unlist(sapply(glist,edgeL), recursive = FALSE, use.names = TRUE)
namesE <- gsub(sprintf("^(%s)\\.",paste(names(glist),collapse="|")),"",names(E))
uniqueE <- lapply(V, function(v) list(edges = unique(unlist(E[namesE==v])) ) )
names(uniqueE) <- V
g <- new("graphNEL", nodes = V, edgeL = uniqueE, edgemode = "directed")
g.bis <- igraph::igraph.from.graphNEL(g)

igraph::plot.igraph(g.bis,vertex.size=2,vertex.label=NA,vertex.color=0,edge.arrow.size=0.2,edge.color=1)

as.list(illuminaHumanv2PATH2PROBE)[gsub("hsa:","",V)]
kegg2pid
unlist(as.list(illuminaHumanv2PMID2PROBE)[gsub("hsa:","",V)])

eid <- unlist(as.list(illuminaHumanv2ENTREZID))
sum(gsub("hsa:","",V) %in% eid)
sub_vec_illumina_id <- eid[sapply( gsub("hsa:","",V), function(id) which(id == eid)[1] )]

# Old annotations : glycolysis / glycogenesis --> wrong networks
f1 <- "/volatile/PACOSE/KEGG/map00030.xml"
f2 <- "/volatile/PACOSE/KEGG/map00230.xml"
f3 <- "/volatile/PACOSE/KEGG/map00010.xml"
# Files downloaded manually at http://rsat.ulb.ac.be/data/KEGG/reference/
g1 <- parseKGML2Graph(f1,genesOnly=FALSE)
g2 <- parseKGML2Graph(f2,genesOnly=FALSE)
g3 <- parseKGML2Graph(f3,genesOnly=FALSE)
glist <- list(map00030=g1,map00230=g2,map00010=g3)

V <- unique(as.vector(unlist(sapply(glist, nodes))))
E <- unlist(sapply(glist,edges), recursive = FALSE, use.names = TRUE)
namesE <- gsub(sprintf("^(%s)\\.",paste(names(glist),collapse="|")),"",names(E))
uniqueE <- lapply(V, function(v) list(edges = unique(unlist(E[namesE==v])) ) )
names(uniqueE) <- V
g <- new("graphNEL", nodes = V, edgeL = uniqueE, edgemode = "directed")
g.bis <- igraph::igraph.from.graphNEL(g)

igraph::plot.igraph(g.bis,vertex.size=2,vertex.label=NA,vertex.color=0,edge.arrow.size=0.2,edge.color=1)

