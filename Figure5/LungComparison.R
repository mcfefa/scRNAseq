install.packages("Seurat")
install.packages("rdetools")

library(Matrix)
library(dplyr)
library(Seurat)
library(readr)
library(rdetools)
library(data.table)

##### Step 1: Load all the datasets
# lung patient 1
lung1norm.data <- Read10X(data.dir = "./data/Patient1/normal/outs/filtered_gene_bc_matrices/GRCh38")
lung1norm <- CreateSeuratObject(raw.data = lung1norm.data, project = "lung1norm")
lung1core.data <- Read10X(data.dir = "./data/Patient1/core/outs/filtered_gene_bc_matrices/GRCh38")
lung1core <- CreateSeuratObject(raw.data = lung1core.data, project = "lung1core")
lung1middle.data <- Read10X(data.dir = "./data/Patient1/middle/outs/filtered_gene_bc_matrices/GRCh38")
lung1middle <- CreateSeuratObject(raw.data = lung1middle.data, project = "lung1middle")
lung1edge.data <- Read10X(data.dir = "./data/Patient1/edge/outs/filtered_gene_bc_matrices/GRCh38")
lung1edge <- CreateSeuratObject(raw.data = lung1edge.data, project = "lung1edge")
# lung patient 2
lung2norm.data <- Read10X(data.dir = "./data/Patient2/normal/outs/filtered_gene_bc_matrices/GRCh38")
lung2norm <- CreateSeuratObject(raw.data = lung2norm.data, project = "lung2norm")
lung2core.data <- Read10X(data.dir = "./data/Patient2/core/outs/filtered_gene_bc_matrices/GRCh38")
lung2core <- CreateSeuratObject(raw.data = lung2core.data, project = "lung2core")
lung2middle.data <- Read10X(data.dir = "./data/Patient2/middle/outs/filtered_gene_bc_matrices/GRCh38")
lung2middle <- CreateSeuratObject(raw.data = lung2middle.data, project = "lung2middle")
lung2edge.data <- Read10X(data.dir = "./data/Patient2/edge/outs/filtered_gene_bc_matrices/GRCh38")
lung2edge <- CreateSeuratObject(raw.data = lung2edge.data, project = "lung2edge")
# lung patient 3
lung3norm.data <- Read10X(data.dir = "./data/Patient3/normal/outs/filtered_gene_bc_matrices/GRCh38")
lung3norm <- CreateSeuratObject(raw.data = lung3norm.data, project = "lung3norm")
lung3core.data <- Read10X(data.dir = "./data/Patient3/core/outs/filtered_gene_bc_matrices/GRCh38")
lung3core <- CreateSeuratObject(raw.data = lung3core.data, project = "lung3core")
lung3middle.data <- Read10X(data.dir = "./data/Patient3/middle/outs/filtered_gene_bc_matrices/GRCh38")
lung3middle <- CreateSeuratObject(raw.data = lung3middle.data, project = "lung3middle")
lung3edge.data <- Read10X(data.dir = "./data/Patient3/edge/outs/filtered_gene_bc_matrices/GRCh38")
lung3edge <- CreateSeuratObject(raw.data = lung3edge.data, project = "lung3edge")
# lung patient 4
lung4norm.data <- Read10X(data.dir = "./data/Patient4/normal/outs/filtered_gene_bc_matrices/GRCh38")
lung4norm <- CreateSeuratObject(raw.data = lung4norm.data, project = "lung4norm")
lung4core.data <- Read10X(data.dir = "./data/Patient4/core/outs/filtered_gene_bc_matrices/GRCh38")
lung4core <- CreateSeuratObject(raw.data = lung4core.data, project = "lung4core")
lung4middle.data <- Read10X(data.dir = "./data/Patient4/middle/outs/filtered_gene_bc_matrices/GRCh38")
lung4middle <- CreateSeuratObject(raw.data = lung4middle.data, project = "lung4middle")
lung4edge.data <- Read10X(data.dir = "./data/Patient4/edge/outs/filtered_gene_bc_matrices/GRCh38")
lung4edge <- CreateSeuratObject(raw.data = lung4edge.data, project = "lung4edge")
# lung patient 5
lung5norm.data <- Read10X(data.dir = "./data/Patient5/normal/outs/filtered_gene_bc_matrices/GRCh38")
lung5norm <- CreateSeuratObject(raw.data = lung5norm.data, project = "lung5norm")
lung5core.data <- Read10X(data.dir = "./data/Patient5/core/outs/filtered_gene_bc_matrices/GRCh38")
lung5core <- CreateSeuratObject(raw.data = lung5core.data, project = "lung5core")
lung5middle.data <- Read10X(data.dir = "./data/Patient5/middle/outs/filtered_gene_bc_matrices/GRCh38")
lung5middle <- CreateSeuratObject(raw.data = lung5middle.data, project = "lung5middle")
lung5edge.data <- Read10X(data.dir = "./data/Patient5/edge/outs/filtered_gene_bc_matrices/GRCh38")
lung5edge <- CreateSeuratObject(raw.data = lung5edge.data, project = "lung5edge")
# lung patient 6
lung6norm.data <- Read10X(data.dir = "./data/Patient6/normal/outs/filtered_gene_bc_matrices/GRCh38")
lung6norm <- CreateSeuratObject(raw.data = lung6norm.data, project = "lung6norm")
lung6core.data <- Read10X(data.dir = "./data/Patient6/core/outs/filtered_gene_bc_matrices/GRCh38")
lung6core <- CreateSeuratObject(raw.data = lung6core.data, project = "lung6core")
lung6edge.data <- Read10X(data.dir = "./data/Patient6/edge/outs/filtered_gene_bc_matrices/GRCh38")
lung6edge <- CreateSeuratObject(raw.data = lung6edge.data, project = "lung6edge")

##### Step 2: Aggregate Datafiles
# aggregate iteratively using MergeSeurat()
lung.combined1 <- MergeSeurat(object1 = lung1norm, object2 = lung1core, add.cell.id1 = "lung1norm", add.cell.id2 = "lung1core", project = "LungComparison")
lung.combined2 <- MergeSeurat(object1 = lung.combined1, object2 = lung1middle, add.cell.id2 = "lung1middle", project = "LungComparison")
lung.combined3 <- MergeSeurat(object1 = lung.combined2, object2 = lung1edge, add.cell.id2 = "lung1edge", project = "LungComparison")
lung.combined4 <- MergeSeurat(object1 = lung.combined3, object2 = lung2norm, add.cell.id2 = "lung2norm", project = "LungComparison")
lung.combined5 <- MergeSeurat(object1 = lung.combined4, object2 = lung2core, add.cell.id2 = "lung2core", project = "LungComparison")
lung.combined6 <- MergeSeurat(object1 = lung.combined5, object2 = lung2middle, add.cell.id2 = "lung2middle", project = "LungComparison")
lung.combined7 <- MergeSeurat(object1 = lung.combined6, object2 = lung2edge, add.cell.id2 = "lung2edge", project = "LungComparison")
lung.combined8 <- MergeSeurat(object1 = lung.combined7, object2 = lung3norm, add.cell.id2 = "lung3norm", project = "LungComparison")
lung.combined9 <- MergeSeurat(object1 = lung.combined8, object2 = lung3core, add.cell.id2 = "lung3core", project = "LungComparison")
lung.combined10 <- MergeSeurat(object1 = lung.combined9, object2 = lung3middle, add.cell.id2 = "lung3middle", project = "LungComparison")
lung.combined11 <- MergeSeurat(object1 = lung.combined10, object2 = lung3edge, add.cell.id2 = "lung3edge", project = "LungComparison")
lung.combined12 <- MergeSeurat(object1 = lung.combined11, object2 = lung4norm, add.cell.id2 = "lung4norm", project = "LungComparison")
lung.combined13 <- MergeSeurat(object1 = lung.combined12, object2 = lung4core, add.cell.id2 = "lung4core", project = "LungComparison")
lung.combined14 <- MergeSeurat(object1 = lung.combined13, object2 = lung4middle, add.cell.id2 = "lung4middle", project = "LungComparison")
lung.combined15 <- MergeSeurat(object1 = lung.combined14, object2 = lung4edge, add.cell.id2 = "lung4edge", project = "LungComparison")
lung.combined16 <- MergeSeurat(object1 = lung.combined15, object2 = lung5norm, add.cell.id2 = "lung5norm", project = "LungComparison")
lung.combined17 <- MergeSeurat(object1 = lung.combined16, object2 = lung5core, add.cell.id2 = "lung5core", project = "LungComparison")
lung.combined18 <- MergeSeurat(object1 = lung.combined17, object2 = lung5middle, add.cell.id2 = "lung5middle", project = "LungComparison")
lung.combined19 <- MergeSeurat(object1 = lung.combined18, object2 = lung5edge, add.cell.id2 = "lung5edge", project = "LungComparison")
lung.combined20 <- MergeSeurat(object1 = lung.combined19, object2 = lung6norm, add.cell.id2 = "lung6norm", project = "LungComparison")
lung.combined21 <- MergeSeurat(object1 = lung.combined20, object2 = lung6core, add.cell.id2 = "lung6core", project = "LungComparison")
lung.combined <- MergeSeurat(object1 = lung.combined21, object2 = lung6edge, add.cell.id2 = "lung6edge", project = "LungComparison")

##### Step 3: Quality Control Step
# identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = lung.combined@data), value = TRUE)
percent.mito <- Matrix::colSums(lung.combined@raw.data[mito.genes, ])/Matrix::colSums(lung.combined@raw.data)
lung.combined <- AddMetaData(object = lung.combined, metadata = percent.mito, col.name = "percent.mito")
lung.combined <- FilterCells(object = lung.combined, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.1))

##### Step 4: Normalize the data
lung.combined <- NormalizeData(object = lung.combined, normalization.method = "LogNormalize", scale.factor = 10000)

##### Step 5: Determine variable genes across cells
lung.combined <- FindVariableGenes(object = lung.combined, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = lung.combined@var.genes)

##### Step 6: Scale data and remove unwanted variation
lung.combined <- ScaleData(object = lung.combined, vars.to.regress = c("nUMI", "percent.mito"), do.par=TRUE, num.cores=2)

##### Step 7: Perform linear dimensional reduction
lung.combined <- RunPCA(object = lung.combined, pc.genes = lung.combined@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
lung.combined <- ProjectPCA(object = lung.combined, do.print = FALSE)

##### Step 8: Cluster Data - Graph-Based Cluster (like Loupe Browser)
# following along: https://satijalab.org/seurat/pbmc3k_tutorial.html 
lung.combinedGraphedBased <- FindClusters(object = lung.combined, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = lung.combinedGraphedBased)
nlevels(lung.combinedGraphedBased@ident)
lung.combinedGraphedBased <- RunTSNE(object = lung.combinedGraphedBased, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = lung.combinedGraphedBased)

###### Step 9: Diversity Scoring
# tested with the graph-based clustering
cat(lung.combinedGraphedBased@cell.names, file="./outputDataNamesLung.csv", sep=",\n")
cat(lung.combinedGraphedBased@ident, file="./outputDataLung.csv", sep=",\n")
### switch to local R studio for plotting
