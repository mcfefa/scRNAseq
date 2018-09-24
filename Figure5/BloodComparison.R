
install.packages("Seurat")
install.package("rdetools")

library(Matrix)
library(dplyr)
library(Seurat)
library(readr)
library(rdetools)
library(data.table)

##### Step 1: Load all the datasets
# Navigate to the directory containing the three files that describe the gene/cell matrix (barcodes.tsv, genes.tsv, matrix.mtx)
CD34.data <- Read10X(data.dir = "./data/CD34/outs/filtered_matrices_mex/hg19/")
CD34 <- CreateSeuratObject(raw.data = CD34.data, project = "CD34+")
CD4.data <- Read10X(data.dir = "./data/CD4/outs/filtered_matrices_mex/hg19/")
CD4 <- CreateSeuratObject(raw.data = CD4.data, project = "CD4+")
CD14.data <- Read10X(data.dir = "./data/CD14/outs/filtered_matrices_mex/hg19/")
CD14 <- CreateSeuratObject(raw.data = CD14.data, project = "CD14+")
CD19.data <- Read10X(data.dir = "./data/CD19/outs/filtered_matrices_mex/hg19/")
CD19 <- CreateSeuratObject(raw.data = CD19.data, project = "CD19+")

##### Step 2: Aggregate Datafiles
# aggregate iteratively using MergeSeurat()
blood.combined1 <- MergeSeurat(object1 = CD34, object2 = CD4, add.cell.id1 = "CD34+", add.cell.id2 = "CD4+", project = "BloodComparison")
blood.combined2 <- MergeSeurat(object1 = blood.combined1, object2 = CD14, add.cell.id2 = "CD14+", project = "BloodComparison")
blood.combined <- MergeSeurat(object1 = blood.combined2, object2 = CD19, add.cell.id2 = "CD19+", project = "BloodComparison")

##### Step 3: Quality Control Step
# identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = blood.combined@data), value = TRUE)
percent.mito <- Matrix::colSums(blood.combined@raw.data[mito.genes, ])/Matrix::colSums(blood.combined@raw.data)
blood.combined <- AddMetaData(object = blood.combined, metadata = percent.mito, col.name = "percent.mito")
blood.combined <- FilterCells(object = blood.combined, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.10))
# cut-offs for nGene identified and Mitochondrial DNA were the same used from: Lambrechts et al, Nature Medicine, 2018

##### Step 4: Normalize the data
blood.combined <- NormalizeData(object = blood.combined, normalization.method = "LogNormalize", scale.factor = 10000)

##### Step 5: Determine variable genes across cells
blood.combined <- FindVariableGenes(object = blood.combined, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = blood.combined@var.genes)

##### Step 6: Scale data and remove unwanted variation
blood.combined <- ScaleData(object = blood.combined, vars.to.regress = c("nUMI", "percent.mito"), do.par=TRUE, num.cores=2)

##### Step 7: Perform linear dimensional reduction
blood.combined <- RunPCA(object = blood.combined, pc.genes = blood.combined@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
blood.combined <- ProjectPCA(object = blood.combined, do.print = FALSE)

##### Step 8: Cluster Data - Graph-Based Cluster (like Loupe Browser)
# adapted from: https://satijalab.org/seurat/pbmc3k_tutorial.html 
blood.combinedGraphedBased <- FindClusters(object = blood.combined, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = blood.combinedGraphedBased)
nlevels(blood.combinedGraphedBased@ident)
blood.combinedGraphedBased <- RunTSNE(object = blood.combinedGraphedBased, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = blood.combinedGraphedBased)

###### Step 9: Diversity Scoring
# collect number of cells per cluster per type
cat(blood.combinedGraphedBased@cell.names, file="./outputDataNamesBlood.csv", sep=",\n")
cat(blood.combinedGraphedBased@ident, file="./outputDataBlood.csv", sep=",\n")
mydat1 <- read.csv("/home/4467528/outputDataBlood.csv")
mydat2 <- read.csv("/home/4467528/outputDataNamesBlood.csv")
fulldat <- cbind(mydat2[1],mydat1[1])
fulltab <- as.data.table(fulldat)
names(fulltab)[1] <- paste("UMI")
names(fulltab)[2] <- paste("cluster")
# group cells by cluster
group_by(fulltab, cluster)
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# then group by type
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)
tabPerClusType <- fulltab %>% group_by(cluster, type) %>% count()
tabPerType <- fulltab %>% group_by(type) %>% count()
# covert the cell counts to frequencies used to calculate the diveristy score
tmp <- matrix(ncol=1, nrow=dim(tabPerClusType)[1])
for(i in 1:dim(tabPerClusType)[1]){
  if(tabPerClusType$type[i]==tabPerType$type[1]) {
    ## if CD34
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[1]
  } else if (tabPerClusType$type[i]==tabPerType$type[2]) {
    ## if CD4
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[2]
  } else if (tabPerClusType$type[i]==tabPerType$type[3]) {
    ## if CD14
     tmp[i] <- tabPerClusType$n[i]/tabPerType$n[3]
  } else if (tabPerClusType$type[i]==tabPerType$type[4]) {
    ## if CD19
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[4]
  } else {
    tmp[i] <- 0.0
  }
}
tmp <- as.data.frame(tmp)
names(tmp)[1] <- paste("freq")
tabPerClusType <- cbind(as.data.frame(tabPerClusType), tmp)

# filter/collect per cell type
CD34clus <- filter(tabPerClusType,type=="CD34+")
CD4clus <- filter(tabPerClusType,type=="CD4+")
CD14clus <- filter(tabPerClusType,type=="CD14+")
CD19clus <- filter(tabPerClusType,type=="CD19+")

# Function to calculate the divesrity index
calcqD <- function(dataset, q){
  diversity <- 0.0;
  for(row in 1:dim(dataset)[1]){
    diversity <- diversity + (dataset$freq[row])^q
  }
  diversity <- diversity^(1/(1-q))
}

# Calculate the diversity index over a range of q (set by qRange)
qRange <- logspace(-2, 2, n = 1000)
qDCD34 <- calcqD(CD34clus,qRange)
qDCD34tab <- as.data.frame(cbind(qRange,qDCD34))
qDCD4 <- calcqD(CD4clus,qRange)
qDCD4tab <- as.data.frame(cbind(qRange,qDCD4))
qDCD14 <- calcqD(CD14clus,qRange)
qDCD14tab <- as.data.frame(cbind(qRange,qDCD14))
qDCD19 <- calcqD(CD19clus,qRange)
qDCD19tab <- as.data.frame(cbind(qRange,qDCD19))

# Plots the spectrum of diversity scores
ggplot(data=qDCD34tab, aes(x=qRange, y=qDCD34, color="CD34")) +
  geom_smooth(size=6) + 
  scale_x_log10() +
  geom_smooth(data=qDCD4tab, aes(x=qRange, y=qDCD4, color="CD4"),size=6) + #, color="navy"
  geom_smooth(data=qDCD14tab, aes(x=qRange, y=qDCD14, color="CD14"),size=6) + #, color="blue"
  geom_smooth(data=qDCD19tab, aes(x=qRange, y=qDCD19, color="CD19"),size=6) + #, color="deepskyblue4"
  labs(title="Diversity Score",x="q", y = "qD") +
  scale_color_manual(name="", 
                     breaks = c("CD14","CD19","CD34","CD4"),
                     values = c("blue","deepskyblue4","dark gray","navy"),
                     guide = guide_legend(override.aes=aes(fill=NA))) +
  theme(rect = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
