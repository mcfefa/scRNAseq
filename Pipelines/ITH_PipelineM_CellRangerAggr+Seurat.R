install.packages("Seurat")
install.packages("rdetools")

library(Matrix)
library(Seurat)
library(dplyr)
library(readr)
library(rdetools)
library(data.table)

##### Step 1: Load all output from cellranger aggr (aggregated AML datasets)
# Navigate to the directory containing the three files that describe the gene/cell matrix (barcodes.tsv, genes.tsv, matrix.mtx)
bmmcAll.data <- Read10X(data.dir = "./data/aggr/outs/filtered_gene_bc_matrices_mex/GRCh38/")
bmmcAll <- CreateSeuratObject(raw.data = bmmcAll.data, project = "AMLvPost")
dim(bmmcAll.data)

##### Step 2: Quality Control Step
# identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = bmmcAll@data), value = TRUE)
percent.mito <- Matrix::colSums(bmmcAll@raw.data[mito.genes, ])/Matrix::colSums(bmmcAll@raw.data)
bmmcAll <- AddMetaData(object = bmmcAll, metadata = percent.mito, col.name = "percent.mito")
bmmcAll <- FilterCells(object = bmmcAll, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.1))
# cut-offs for nGene identified and Mitochondrial DNA were the same used from: Lambrechts et al, Nature Medicine, 2018

##### Step 4: Normalize the data
bmmcAll <- NormalizeData(object = bmmcAll, normalization.method = "LogNormalize", scale.factor = 10000)

##### Step 5: Determine variable genes across cells
bmmcAll <- FindVariableGenes(object = bmmcAll, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = bmmcAll@var.genes)

##### Step 6: Scale data and remove unwanted variation
bmmcAll <- ScaleData(object = bmmcAll, vars.to.regress = c("nUMI", "percent.mito"), do.par=TRUE, num.cores=2)

##### Step 7: Perform linear dimensional reduction
bmmcAll <- RunPCA(object = bmmcAll, pc.genes = bmmcAll@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = bmmcAll, dim.1 = 1, dim.2 = 2)
bmmcAll <- ProjectPCA(object = bmmcAll, do.print = FALSE)

##### Step 8: Cluster Data - Graph-Based Cluster (like Loupe Browser)
# following along: https://satijalab.org/seurat/pbmc3k_tutorial.html 
bmmcAll.combinedGraphedBased <- FindClusters(object = bmmcAll, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = bmmcAll.combinedGraphedBased)
nlevels(bmmcAll.combinedGraphedBased@ident)
bmmcAll.combinedGraphedBased <- RunTSNE(object = bmmcAll.combinedGraphedBased, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = bmmcAll.combinedGraphedBased)
dev.off()

###### Step 9: Diversity Scoring
# collect number of cells per cluster per type
cat(bmmcAll.combinedGraphedBased@cell.names, file="./outputDataNamesAML.csv", sep=",\n")
cat(bmmcAll.combinedGraphedBased@ident, file="./outputDataAML.csv", sep=",\n")
mydat1 <- read.csv("./outputDataAML.csv")
mydat2 <- read.csv("./outputDataNamesAML.csv")
fulldat <- cbind(mydat2[1],mydat1[1])
# count number of unique barcodes per cluster
fulltab <- as.data.table(fulldat)
names(fulltab)[1] <- paste("UMI")
names(fulltab)[2] <- paste("cluster")
# group cells by cluster
group_by(fulltab, cluster)
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# then group by type
type <- sub(".*-","",fulltab$UMI)
fulltab <- cbind(fulltab, type)
tabPerClusType <- fulltab %>% group_by(cluster, type) %>% count()
tabPerType <- fulltab %>% group_by(type) %>% count()
# covert the cell counts to frequencies used to calculate the diveristy score
tmp <- matrix(ncol=1, nrow=dim(tabPerClusType)[1])
for(i in 1:dim(tabPerClusType)[1]){
  if(tabPerClusType$type[i]==tabPerType$type[1]) {
    ## if AML027pre
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[1]
  } else if (tabPerClusType$type[i]==tabPerType$type[2]) {
    ## if AML027post
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[2]
  } else if (tabPerClusType$type[i]==tabPerType$type[3]) {
    ## if AML035pre
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[3]
  } else if (tabPerClusType$type[i]==tabPerType$type[4]) {
    ## if AML035post
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[4]
  } else if (tabPerClusType$type[i]==tabPerType$type[5]) {
    ## if Healthy 1
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[5]
  } else if (tabPerClusType$type[i]==tabPerType$type[6]) {
    ## if Healthy 2
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[6]
  }else {
    tmp[i] <- 0.0
  }
}
tmp <- as.data.frame(tmp)
names(tmp)[1] <- paste("freq")
tabPerClusType <- cbind(as.data.frame(tabPerClusType), tmp)

# filter/collect per cell type
AML1 <- filter(tabPerClusType,type==1)#"027pre")
AML2 <- filter(tabPerClusType,type==3)#"035pre")
H1 <- filter(tabPerClusType,type==5)#"H1")
H2 <- filter(tabPerClusType,type==6)#"H2")
postAML1 <- filter(tabPerClusType,type==2)#"027post")
postAML2 <- filter(tabPerClusType,type==4)#"035post")

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
qDAML1 <- calcqD(AML1,qRange)
qDAML1tab <- as.data.frame(cbind(qRange,qDAML1))
qDAML2 <- calcqD(AML2,qRange)
qDAML2tab <- as.data.frame(cbind(qRange,qDAML2))

qDH1 <- calcqD(H1,qRange)
qDH1tab <- as.data.frame(cbind(qRange,qDH1))
qDH2 <- calcqD(H2,qRange)
qDH2tab <- as.data.frame(cbind(qRange,qDH2))

qDpostAML1 <- calcqD(postAML1,qRange)
qDpostAML1tab <- as.data.frame(cbind(qRange,qDpostAML1))
qDpostAML2 <- calcqD(postAML2,qRange)
qDpostAML2tab <- as.data.frame(cbind(qRange,qDpostAML2))

### Plots the spectrum of diversity scores
# Plotting all data 
ggplot(data=qDAML1tab, aes(x=qRange, y=qDAML1)) +
  geom_smooth(color="dark red") + 
  scale_x_log10() +
  geom_smooth(data=qDAML2tab, aes(x=qRange, y=qDAML2), color="red") +
  geom_smooth(data=qDH1tab, aes(x=qRange, y=qDH1), color="dark gray") + 
  geom_smooth(data=qDH2tab, aes(x=qRange, y=qDH2), color="black") +
  labs(title="Diversity Score",x="q", y = "qD") + 
  geom_smooth(data=qDpostAML2tab, aes(x=qRange, y=qDpostAML1), color="purple") + 
  geom_smooth(data=qDpostAML2tab, aes(x=qRange, y=qDpostAML2), color="blueviolet")

# Plotting AMLpre vs AMLpost
ggplot(data=qDAML1tab, aes(x=qRange, y=qDAML1)) +
  geom_smooth(color="dark red") + 
  scale_x_log10() +
  geom_smooth(data=qDAML2tab, aes(x=qRange, y=qDAML2), color="red") +
  labs(title="Diversity Score",x="q", y = "qD") + 
  geom_smooth(data=qDpostAML2tab, aes(x=qRange, y=qDpostAML1), color="purple") + 
  geom_smooth(data=qDpostAML2tab, aes(x=qRange, y=qDpostAML2), color="blueviolet")

# Plotting AMLpre vs Healthy
ggplot(data=qDAML1tab, aes(x=qRange, y=qDAML1)) +
  geom_smooth(color="dark red") + 
  scale_x_log10() +
  geom_smooth(data=qDAML2tab, aes(x=qRange, y=qDAML2), color="red") +
  geom_smooth(data=qDH1tab, aes(x=qRange, y=qDH1), color="dark gray") + 
  geom_smooth(data=qDH2tab, aes(x=qRange, y=qDH2), color="black") +
  labs(title="Diversity Score",x="q", y = "qD") 

