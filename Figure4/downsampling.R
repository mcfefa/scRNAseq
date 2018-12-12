######## 9/2018
########  Automating Downsampling of 10x scRNAseq data
#           Step 1: Prior to downsampling, create Seurat object and save
#           Step 2: Load Seurat object
#           Step 3: Subsampling function
#           Step 4: Record clusters found after subsampling
#           Step 5: Calculate and record stats on clusters found from subset

## Load necessary libraries
library(Matrix)
library(readr)
library(Seurat)
library(dplyr)
library(readr)
library(rdetools)
library(data.table)

## Set initial inputs/arguments
dataFrom <- # add string of cell type in double quotes, ex: "AML027pre"
datasetFile <- # add path to rds file for cell type subsetting, ex: "./data/AML027pre_20180920.rds"
per <- # add here the fraction of cells you'd like to subset, ex: 0.10 indiciates a subsample to 10% the total number of cells in the experiment
runs <- # add the numbre of times you'd like to repeat this, ex: 1000
fileName <- paste(dataFrom,"_per",per*100,"runs",runs,sep="")
tmpfile1 <- paste("./outputData_per",per*100,".csv",sep="");
tmpfile2 <- paste("./outputDataNames_per",per*100,".csv",sep="");

## Load Seurat dataset object
bmmc <- readRDS(datasetFile)

## Define function to subset and cluster Seurat dataset object
#    Inputs: dataset = Seurat object (quality controlled, normalized 10x dataset)
#            sampPer = fraction of the number of cells interested in sampling (ie 10% of 100 cells will return number of clusters from a 10-cell subset)
#            ranNum = based on iteration number, set to make sure don't continously produce the same subsample at the same sampPer cutoff
#
subClust <- function(dataset,sampPer,ranNum) {
  numCells <- dim(dataset@data)[2];
  subset <- SubsetData(dataset, cells.use = sample(x = dataset@cell.names, size = numCells*sampPer),random.seed = ranNum );
  subset <- ScaleData(object = subset, vars.to.regress = c("nUMI","percent.mito"), display.progress=FALSE);
  subset <- RunPCA(object = subset, pc.genes = subset@var.genes, do.print = FALSE);
  subset <- ProjectPCA(object = subset, do.print = FALSE);
  subset <- FindClusters(object = subset, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc=TRUE);
  #PrintFindClustersParams(object = subset);
  subset <- RunTSNE(object = subset, dims.use = 1:20, do.fast = TRUE);
  # array of levels: levels(subset@ident)
  clus <- nlevels(subset@ident)
  
  cat(subset@cell.names, file=tmpfile2, sep=",\n", append=FALSE)
  cat(subset@ident, file=tmpfile1, sep=",\n", append=FALSE)
  mydat1 <- read.csv(tmpfile1)
  mydat2 <- read.csv(tmpfile2)
  
  fulldat <- cbind(mydat2[1],mydat1[1])
  fulltab <- as.data.table(fulldat)
  names(fulltab)[1] <- paste("UMI")
  names(fulltab)[2] <- paste("cluster")
  group_by(fulltab, cluster)
  tabPerClus <- fulltab %>% group_by(cluster) %>% count()
  
  tmp <- matrix(ncol=1, nrow=dim(tabPerClus)[1])
  for(i in 1:dim(tabPerClus)[1]){
      tmp[i] <- tabPerClus$n[i]/dim(subset@data)[2]
  }
  tmp <- as.data.frame(tmp)
  names(tmp)[1] <- paste("freq")
  tabPerClusFreq <- cbind(as.data.frame(tabPerClus), tmp)

  calcqD <- function(dat, q){
    diversity <- 0.0;
    if(q == 1.0){
      for(row in 1:dim(dat)[1]){
          diversity <- diversity + log((dat$freq[row])^(dat$freq[row]))
          }
      diversity <- abs(diversity)
    }else{
      for(row in 1:dim(dat)[1]){
        diversity <- diversity + (dat$freq[row])^q
      }
      diversity <- diversity^(1/(1-q))
    }
  }
  
  qDp01 <- calcqD(tabPerClusFreq,0.01)
  qDp1 <- calcqD(tabPerClusFreq,0.1)
  qD1 <- calcqD(tabPerClusFreq,1.0)
  qD2 <- calcqD(tabPerClusFreq,2.0)
  qD10 <- calcqD(tabPerClusFreq,10.0)
  qD100 <- calcqD(tabPerClusFreq,100.0)
  
  cbind(clus,qDp01,qDp1,qD1,qD2,qD10,qD100)
}

## Set up subsampling based on initial conditions previous set (to repeat RUNS number of times)
# per = 0.10, runs = 100 would subsample the dataset to 10% of cells 100 different times and 
#                        record the number of clusters found in each 10% subsampling run
output <- matrix(ncol=7, nrow=runs)
for (i in 1:runs){
  output[i,] <- subClust(bmmc,per,i)
  print(i)
}
# write output to a file
write.csv(output, paste(fileName,".csv",sep=""))

