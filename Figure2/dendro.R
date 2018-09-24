
install.packages("tripack")
install.packages("colorRamps")
install.packages("qgraph")
install.packages("plotrix")


library("colorRamps")
library(tripack)
library(RColorBrewer)



library(qgraph)
library("plotrix")


mat <- read.csv("DistanceMatrix_MeanUMI.csv",header=T)
graph <- hclust(as.dist(mat),  method = "cen")
plot(graph,hang = -1, cex = 0.8, ylab = "distance")

dist_mi <- 1/as.matrix(mat)

size <- as.matrix(read.csv("clustersizes.csv",header=F))

colscale <- as.matrix(read.csv("colorListHEX.csv",header=F))

png("Fig2D_dendrogram_out.png",width=16,height=12,units = "cm",res=1200)
qgraph(dist_mi, layout='spring', vsize=size , color=colscale , label.cex=1.2) 
dev.off()

