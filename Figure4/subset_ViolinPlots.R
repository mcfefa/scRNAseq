######## 9/2018
########  Automating Downsampling of 10x scRNAseq data
#           Step 1: Load CSV files of subsampled data
#           Step 2: Create a dataframe across all subsampled fractions
#           Step 3: Plot distribution of qD scores 

# load libraries
library(readr)
library(ggplot2)

## load data in 
# after subsampling, read in each csv with subsamples results like shown below 
# previously subsampled files are included in the github Figure4/data/Subset-CSV-files 
fileName <- "./data/Healthy1/Healthy1_per50runs1000"
Healthy1per50runs1k <- read_csv(paste(fileName,".csv",sep=""))
names(Healthy1per50runs1k)[2] <- paste("NumClus")
names(Healthy1per50runs1k)[3] <- paste("q=0.01")
names(Healthy1per50runs1k)[4] <- paste("q=0.1")
names(Healthy1per50runs1k)[5] <- paste("q=1")
names(Healthy1per50runs1k)[6] <- paste("q=2")
names(Healthy1per50runs1k)[7] <- paste("q=10")
names(Healthy1per50runs1k)[8] <- paste("q=100")

fileName <- "./data/Healthy2/Healthy2_per50runs1000"
Healthy2per50runs1k <- read_csv(paste(fileName,".csv",sep=""))
names(Healthy2per50runs1k)[2] <- paste("NumClus")
names(Healthy2per50runs1k)[3] <- paste("q=0.01")
names(Healthy2per50runs1k)[4] <- paste("q=0.1")
names(Healthy2per50runs1k)[5] <- paste("q=1")
names(Healthy2per50runs1k)[6] <- paste("q=2")
names(Healthy2per50runs1k)[7] <- paste("q=10")
names(Healthy2per50runs1k)[8] <- paste("q=100")

fileName <- "./data/AML027pre/AML027pre_per50runs1000"
preAML027per50runs1k <- read_csv(paste(fileName,".csv",sep=""))
names(preAML027per50runs1k)[2] <- paste("NumClus")
names(preAML027per50runs1k)[3] <- paste("q=0.01")
names(preAML027per50runs1k)[4] <- paste("q=0.1")
names(preAML027per50runs1k)[5] <- paste("q=1")
names(preAML027per50runs1k)[6] <- paste("q=2")
names(preAML027per50runs1k)[7] <- paste("q=10")
names(preAML027per50runs1k)[8] <- paste("q=100")

fileName <- "./data/AML027post/AML027post_per50runs1000"
postAML027per50runs1k <- read_csv(paste(fileName,".csv",sep=""))
names(postAML027per50runs1k)[2] <- paste("NumClus")
names(postAML027per50runs1k)[3] <- paste("q=0.01")
names(postAML027per50runs1k)[4] <- paste("q=0.1")
names(postAML027per50runs1k)[5] <- paste("q=1")
names(postAML027per50runs1k)[6] <- paste("q=2")
names(postAML027per50runs1k)[7] <- paste("q=10")
names(postAML027per50runs1k)[8] <- paste("q=100")

fileName <- "./data/AML035pre/AML035pre_per50runs1000"
preAML035per50runs1k <- read_csv(paste(fileName,".csv",sep=""))
names(preAML035per50runs1k)[2] <- paste("NumClus")
names(preAML035per50runs1k)[3] <- paste("q=0.01")
names(preAML035per50runs1k)[4] <- paste("q=0.1")
names(preAML035per50runs1k)[5] <- paste("q=1")
names(preAML035per50runs1k)[6] <- paste("q=2")
names(preAML035per50runs1k)[7] <- paste("q=10")
names(preAML035per50runs1k)[8] <- paste("q=100")

fileName <- "./data/AML035post/AML035post_per50runs1000"
postAML035per50runs1k <- read_csv(paste(fileName,".csv",sep=""))
names(postAML035per50runs1k)[2] <- paste("NumClus")
names(postAML035per50runs1k)[3] <- paste("q=0.01")
names(postAML035per50runs1k)[4] <- paste("q=0.1")
names(postAML035per50runs1k)[5] <- paste("q=1")
names(postAML035per50runs1k)[6] <- paste("q=2")
names(postAML035per50runs1k)[7] <- paste("q=10")
names(postAML035per50runs1k)[8] <- paste("q=100")

dataset <- read_csv(paste(fileName,".csv",sep=""))
names(dataset)[2] <- paste("NumClus")
names(dataset)[3] <- paste("q=0.01")
names(dataset)[4] <- paste("q=0.1")
names(dataset)[5] <- paste("q=1")
names(dataset)[6] <- paste("q=2")
names(dataset)[7] <- paste("q=10")
names(dataset)[8] <- paste("q=100")

# create metrics file
metricFileName <- paste(fileName,"_qD_metrics.txt",sep="")
cat("Fraction Subsampled: ", file=metricFileName, sep="")
cat(fileName, file=metricFileName, sep="\n", append=TRUE)
cat("Mean number of clusters: ", file=metricFileName, sep="", append=TRUE)
runsMean <- mean(as.matrix(dataset[,2]))
cat(runsMean, file=metricFileName, sep="\n", append=TRUE)
cat("Median number of clusters: ", file=metricFileName, sep="", append=TRUE)
runsMed <- median(as.matrix(dataset[,2]))
cat(runsMed, file=metricFileName, sep="\n", append=TRUE)
cat("Standard deviation of number of clusters: ", file=metricFileName, sep="", append=TRUE)
runsSD <- sd(as.matrix(dataset[,2]))
cat(runsSD, file=metricFileName, sep="\n", append=TRUE)
cat("Minimum number of clusters: ", file=metricFileName, sep="", append=TRUE)
runsMin <- min(as.matrix(dataset[,2]))
cat(runsMin, file=metricFileName, sep="\n", append=TRUE)
cat("Maximum number of clusters: ", file=metricFileName, sep="", append=TRUE)
runsMax <- max(as.matrix(dataset[,2]))
cat(runsMax, file=metricFileName, sep="\n", append=TRUE)

## PLOTTING

# SUBTRACT MEAN FROM EACH SO PLOTTING RELATIVE CHANGES IN qD

# Healthy1-qD-distribution-2018-09-21
#Healthy1_ViolinPlot_1k_relqD_2018-09-22
qD1transH1 <- exp(Healthy1per50runs1k[,5])
qD1H1backup <- Healthy1per50runs1k[,5]
Healthy1per50runs1k[,5] <- qD1transH1
dfhealthy1 <- as.data.frame(Healthy1per50runs1k[,3:8])
Healthy1per50runs1kM1 <- mean(as.matrix(dfhealthy1[,1]))
Healthy1per50runs1kM2 <- mean(as.matrix(dfhealthy1[,2]))
Healthy1per50runs1kM3 <- mean(as.matrix(dfhealthy1[,3]))
Healthy1per50runs1kM4 <- mean(as.matrix(dfhealthy1[,4]))
Healthy1per50runs1kM5 <- mean(as.matrix(dfhealthy1[,5]))
Healthy1per50runs1kM6 <- mean(as.matrix(dfhealthy1[,6]))
dfhealthy1[,1] <- dfhealthy1[,1] - Healthy1per50runs1kM1
dfhealthy1[,2] <- dfhealthy1[,2] - Healthy1per50runs1kM2
dfhealthy1[,3] <- dfhealthy1[,3] - Healthy1per50runs1kM3
dfhealthy1[,4] <- dfhealthy1[,4] - Healthy1per50runs1kM4
dfhealthy1[,5] <- dfhealthy1[,5] - Healthy1per50runs1kM5
dfhealthy1[,6] <- dfhealthy1[,6] - Healthy1per50runs1kM6
dfhealthy1.m <- reshape2::melt(dfhealthy1,id.vars=NULL)
ggplot(dfhealthy1.m, aes(x=variable, y=value, fill=variable)) +
    geom_violin(scale="width",adjust = 1,width = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="white") +
    # scale_fill_brewer(palette="gray38")+
    labs(title="qD Healthy 1",x="", y = "qD") +
    ylim(-4,4) +
    theme_bw()

#Healthy2_ViolinPlot_1k_relqD_2018-09-22
qD1transH2 <- exp(Healthy2per50runs1k[,5])
qD1H2backup <- Healthy2per50runs1k[,5]
Healthy2per50runs1k[,5] <- qD1transH2
dfhealthy2 <- as.data.frame(Healthy2per50runs1k[,3:8])
Healthy2per50runs1kM1 <- mean(as.matrix(dfhealthy2[,1]))
Healthy2per50runs1kM2 <- mean(as.matrix(dfhealthy2[,2]))
Healthy2per50runs1kM3 <- mean(as.matrix(dfhealthy2[,3]))
Healthy2per50runs1kM4 <- mean(as.matrix(dfhealthy2[,4]))
Healthy2per50runs1kM5 <- mean(as.matrix(dfhealthy2[,5]))
Healthy2per50runs1kM6 <- mean(as.matrix(dfhealthy2[,6]))
dfhealthy2[,1] <- dfhealthy2[,1] - Healthy2per50runs1kM1
dfhealthy2[,2] <- dfhealthy2[,2] - Healthy2per50runs1kM2
dfhealthy2[,3] <- dfhealthy2[,3] - Healthy2per50runs1kM3
dfhealthy2[,4] <- dfhealthy2[,4] - Healthy2per50runs1kM4
dfhealthy2[,5] <- dfhealthy2[,5] - Healthy2per50runs1kM5
dfhealthy2[,6] <- dfhealthy2[,6] - Healthy2per50runs1kM6
dfhealthy2.m <- reshape2::melt(dfhealthy2,id.vars=NULL)
ggplot(dfhealthy2.m, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale="width",adjust = 1,width = 0.5) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="white") +
  # scale_fill_brewer(palette="gray38")+
  labs(title="qD Healthy 2",x="", y = "qD") +
  ylim(-4,4) +
  theme_bw()

# AML027pre-qD-distribution-2018-09-21
# AML027pre_ViolinPlot_1k_relqD_2018-09-22
qD1transAML027 <- exp(preAML027per50runs1k[,5])
qD1AML027backup <- preAML027per50runs1k[,5]
preAML027per50runs1k[,5] <- qD1transAML027
dfaml027pre <- as.data.frame(preAML027per50runs1k[,3:8])
preAML027per50runs1kM1 <- mean(as.matrix(dfaml027pre[,1]))
preAML027per50runs1kM2 <- mean(as.matrix(dfaml027pre[,2]))
preAML027per50runs1kM3 <- mean(as.matrix(dfaml027pre[,3]))
preAML027per50runs1kM4 <- mean(as.matrix(dfaml027pre[,4]))
preAML027per50runs1kM5 <- mean(as.matrix(dfaml027pre[,5]))
preAML027per50runs1kM6 <- mean(as.matrix(dfaml027pre[,6]))
dfaml027pre[,1] <- dfaml027pre[,1] - preAML027per50runs1kM1
dfaml027pre[,2] <- dfaml027pre[,2] - preAML027per50runs1kM2
dfaml027pre[,3] <- dfaml027pre[,3] - preAML027per50runs1kM3
dfaml027pre[,4] <- dfaml027pre[,4] - preAML027per50runs1kM4
dfaml027pre[,5] <- dfaml027pre[,5] - preAML027per50runs1kM5
dfaml027pre[,6] <- dfaml027pre[,6] - preAML027per50runs1kM6
dfaml027pre.m <- reshape2::melt(dfaml027pre,id.vars=NULL)
ggplot(dfaml027pre.m, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale="width",adjust = 1,width = 0.5) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="white") +
  # scale_fill_brewer(palette="gray38")+
  labs(title="qD AML027 pre",x="", y = "qD") +
  ylim(-4,4) +
  theme_bw()

# AML035pre-qD-distribution-2018-09-21
# AML035pre_ViolinPlot_1k_relqD_2018-09-22
qD1transAML035 <- exp(preAML035per50runs1k[,5])
qD1AML035backup <- preAML035per50runs1k[,5]
preAML035per50runs1k[,5] <- qD1transAML035
dfaml035pre <- as.data.frame(preAML035per50runs1k[,3:8])
dfaml035pre <- as.data.frame(preAML035per50runs1k[,3:8])
preAML035per50runs1kM1 <- mean(as.matrix(dfaml035pre[,1]))
preAML035per50runs1kM2 <- mean(as.matrix(dfaml035pre[,2]))
preAML035per50runs1kM3 <- mean(as.matrix(dfaml035pre[,3]))
preAML035per50runs1kM4 <- mean(as.matrix(dfaml035pre[,4]))
preAML035per50runs1kM5 <- mean(as.matrix(dfaml035pre[,5]))
preAML035per50runs1kM6 <- mean(as.matrix(dfaml035pre[,6]))
dfaml035pre[,1] <- dfaml035pre[,1] - preAML035per50runs1kM1
dfaml035pre[,2] <- dfaml035pre[,2] - preAML035per50runs1kM2
dfaml035pre[,3] <- dfaml035pre[,3] - preAML035per50runs1kM3
dfaml035pre[,4] <- dfaml035pre[,4] - preAML035per50runs1kM4
dfaml035pre[,5] <- dfaml035pre[,5] - preAML035per50runs1kM5
dfaml035pre[,6] <- dfaml035pre[,6] - preAML035per50runs1kM6
dfaml035pre.m <- reshape2::melt(dfaml035pre,id.vars=NULL)
ggplot(dfaml035pre.m, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale="width",adjust = 1,width = 0.5) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="white") +
  # scale_fill_brewer(palette="gray38")+
  labs(title="qD AML035 pre",x="", y = "qD") +
  ylim(-4,4) +
  theme_bw()

# AML027post-qD-distribution-2018-09-21
# AML027post_ViolinPlot_1k_relqD_2018-09-22
qD1transAML027post <- exp(postAML027per50runs1k[,5])
qD1AML027postbackup <- postAML027per50runs1k[,5]
postAML027per50runs1k[,5] <- qD1transAML027post
dfaml027post <- as.data.frame(postAML027per50runs1k[,3:8])
postAML027per50runs1kM1 <- mean(as.matrix(dfaml027post[,1]))
postAML027per50runs1kM2 <- mean(as.matrix(dfaml027post[,2]))
postAML027per50runs1kM3 <- mean(as.matrix(dfaml027post[,3]))
postAML027per50runs1kM4 <- mean(as.matrix(dfaml027post[,4]))
postAML027per50runs1kM5 <- mean(as.matrix(dfaml027post[,5]))
postAML027per50runs1kM6 <- mean(as.matrix(dfaml027post[,6]))
dfaml027post[,1] <- dfaml027post[,1] - postAML027per50runs1kM1
dfaml027post[,2] <- dfaml027post[,2] - postAML027per50runs1kM2
dfaml027post[,3] <- dfaml027post[,3] - postAML027per50runs1kM3
dfaml027post[,4] <- dfaml027post[,4] - postAML027per50runs1kM4
dfaml027post[,5] <- dfaml027post[,5] - postAML027per50runs1kM5
dfaml027post[,6] <- dfaml027post[,6] - postAML027per50runs1kM6
dfaml027post.m <- reshape2::melt(dfaml027post,id.vars=NULL)
ggplot(dfaml027post.m, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale="width",adjust = 1,width = 0.5) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="white") +
  # scale_fill_brewer(palette="gray38")+
  labs(title="qD AML027 post",x="", y = "qD") +
  ylim(-4,4) +
  theme_bw()

# AML035post-qD-distribution-2018-09-21
# AML035post_ViolinPlot_1k_relqD_2018-09-22
qD1transAML035post <- exp(postAML035per50runs1k[,5])
qD1AML035postbackup <- postAML035per50runs1k[,5]
postAML035per50runs1k[,5] <- qD1transAML035post
dfaml035post <- as.data.frame(postAML035per50runs1k[,3:8])
postAML035per50runs1kM1 <- mean(as.matrix(dfaml035post[,1]))
postAML035per50runs1kM2 <- mean(as.matrix(dfaml035post[,2]))
postAML035per50runs1kM3 <- mean(as.matrix(dfaml035post[,3]))
postAML035per50runs1kM4 <- mean(as.matrix(dfaml035post[,4]))
postAML035per50runs1kM5 <- mean(as.matrix(dfaml035post[,5]))
postAML035per50runs1kM6 <- mean(as.matrix(dfaml035post[,6]))
dfaml035post[,1] <- dfaml035post[,1] - postAML035per50runs1kM1
dfaml035post[,2] <- dfaml035post[,2] - postAML035per50runs1kM2
dfaml035post[,3] <- dfaml035post[,3] - postAML035per50runs1kM3
dfaml035post[,4] <- dfaml035post[,4] - postAML035per50runs1kM4
dfaml035post[,5] <- dfaml035post[,5] - postAML035per50runs1kM5
dfaml035post[,6] <- dfaml035post[,6] - postAML035per50runs1kM6
dfaml035post.m <- reshape2::melt(dfaml035post,id.vars=NULL)
ggplot(dfaml035post.m, aes(x=variable, y=value, fill=variable)) +
  geom_violin(scale="width",adjust = 1,width = 0.5) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill="white") +
  # scale_fill_brewer(palette="gray38")+
  labs(title="qD AML035 post",x="", y = "qD") +
  ylim(-4,4) +
  theme_bw()

