install.packages("Seurat") 
install.packages("rdetools")

library(Matrix)
library(dplyr)
library(Seurat)
library(readr)
library(rdetools)
library(data.table)

###### Step 9: Diversity Scoring
# load data exported by LungComparison.R
mydat1 <- read.csv("./outputDataLung.csv")
mydat2 <- read.csv("./outputDataNamesLung.csv")
fulldat <- cbind(mydat2[1],mydat1[1])
fulltab <- as.data.table(fulldat)
# name table columns
names(fulltab)[1] <- paste("UMI")
names(fulltab)[2] <- paste("cluster")
# group data based on clusters
group_by(fulltab, cluster)
# create a table counting unqiue UMIs/cells per cluster
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)
# creates a table containing the number of cells per cluster per condition
tabPerClusType <- fulltab %>% group_by(cluster, type) %>% count()
# create a table containing the number of unique per condition (lung1norm, etc.)
tabPerType <- fulltab %>% group_by(type) %>% count()

# calculate the frequency of each cells out of total number of cells per condition in each cluster
# these are depended on total number of cells in the specific condition (tabPerType)
# the fraction of cells in cluster X out of all cells in condition Y
tmp <- matrix(ncol=1, nrow=dim(tabPerClusType)[1])
for(i in 1:dim(tabPerClusType)[1]){
  if(tabPerClusType$type[i]==tabPerType$type[1]) {
    ## if lung1norm
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[1]
  } else if (tabPerClusType$type[i]==tabPerType$type[2]) {
    ## if lung1core
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[2]
  } else if (tabPerClusType$type[i]==tabPerType$type[3]) {
    ## if lung1middle
     tmp[i] <- tabPerClusType$n[i]/tabPerType$n[3]
  } else if (tabPerClusType$type[i]==tabPerType$type[4]) {
    ## if lung1edge
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[4]
  } else if (tabPerClusType$type[i]==tabPerType$type[5]) {
    ## if lung2norm
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[5]
  } else if (tabPerClusType$type[i]==tabPerType$type[6]) {
    ## if lung2core
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[6]
  } else if (tabPerClusType$type[i]==tabPerType$type[7]) {
    ## if lung2middle
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[7]
  } else if (tabPerClusType$type[i]==tabPerType$type[8]) {
    ## if lung2edge
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[8]
  } else if (tabPerClusType$type[i]==tabPerType$type[9]) {
    ## if lung3norm
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[9]
  } else if (tabPerClusType$type[i]==tabPerType$type[10]) {
    ## if lung3core
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[10]
  } else if (tabPerClusType$type[i]==tabPerType$type[11]) {
    ## if lung3middle
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[11]
  } else if (tabPerClusType$type[i]==tabPerType$type[12]) {
    ## if lung3edge
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[12]  
  } else if (tabPerClusType$type[i]==tabPerType$type[13]) {
    ## if lung4norm
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[13]
  } else if (tabPerClusType$type[i]==tabPerType$type[14]) {
    ## if lung4core
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[14]
  } else if (tabPerClusType$type[i]==tabPerType$type[15]) {
    ## if lung4middle
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[15]
  } else if (tabPerClusType$type[i]==tabPerType$type[16]) {
    ## if lung4edge
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[16]
  } else if (tabPerClusType$type[i]==tabPerType$type[17]) {
    ## if lung5norm
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[17]
  } else if (tabPerClusType$type[i]==tabPerType$type[18]) {
    ## if lung5core
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[18]
  } else if (tabPerClusType$type[i]==tabPerType$type[19]) {
    ## if lung5middle
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[19]
  } else if (tabPerClusType$type[i]==tabPerType$type[20]) {
    ## if lung5edge
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[20]    
  } else if (tabPerClusType$type[i]==tabPerType$type[21]) {
    ## if lung6norm
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[21]
  } else if (tabPerClusType$type[i]==tabPerType$type[22]) {
    ## if lung6core
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[22]
  } else if (tabPerClusType$type[i]==tabPerType$type[23]) {
    ## if lung6edge
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[23]  
  } else {
    tmp[i] <- 0.0
  }
}
tmp <- as.data.frame(tmp)
names(tmp)[1] <- paste("freq")
tabPerClusType <- cbind(as.data.frame(tabPerClusType), tmp)

# filtering the overall table to have a smaller table with just 1 condition in it
lung1normalClus <- filter(tabPerClusType,type=="lung1norm")
lung1coreClus <- filter(tabPerClusType,type=="lung1core")
lung1middleClus <- filter(tabPerClusType,type=="lung1middle")
lung1edgeClus <- filter(tabPerClusType,type=="lung1edge")
lung2normalClus <- filter(tabPerClusType,type=="lung2norm")
lung2coreClus <- filter(tabPerClusType,type=="lung2core")
lung2middleClus <- filter(tabPerClusType,type=="lung2middle")
lung2edgeClus <- filter(tabPerClusType,type=="lung2edge")
lung3normalClus <- filter(tabPerClusType,type=="lung3norm")
lung3coreClus <- filter(tabPerClusType,type=="lung3core")
lung3middleClus <- filter(tabPerClusType,type=="lung3middle")
lung3edgeClus <- filter(tabPerClusType,type=="lung3edge")
lung4normalClus <- filter(tabPerClusType,type=="lung4norm")
lung4coreClus <- filter(tabPerClusType,type=="lung4core")
lung4middleClus <- filter(tabPerClusType,type=="lung4middle")
lung4edgeClus <- filter(tabPerClusType,type=="lung4edge")
lung5normalClus <- filter(tabPerClusType,type=="lung5norm")
lung5coreClus <- filter(tabPerClusType,type=="lung5core")
lung5middleClus <- filter(tabPerClusType,type=="lung5middle")
lung5edgeClus <- filter(tabPerClusType,type=="lung5edge")
lung6normalClus <- filter(tabPerClusType,type=="lung6norm")
lung6coreClus <- filter(tabPerClusType,type=="lung6core")
lung6edgeClus <- filter(tabPerClusType,type=="lung6edge")

# aggregating individual condition tables based on similar contions
lungnormalClus <- rbind(lung1normalClus,lung2normalClus,lung3normalClus,lung4normalClus,lung5normalClus,lung6normalClus)
lungnormalClus <- select(lungnormalClus, one_of(c("cluster", "n")))
tmpNorm <- matrix(ncol=2, nrow=dim(tabPerClus)[1])
for(j in 1:dim(tabPerClus)[1]){
  iter <- as.numeric(j)
  tmpNorm[j,1] <- iter
  a <- lapply(filter(lungnormalClus, lungnormalClus$cluster == j) %>% summarise_at(c("n"),sum,na.rm=TRUE), as.numeric)
  b <- a[[1]]
  c <- (b/(tabPerType$n[4]+tabPerType$n[8]+tabPerType$n[12]+tabPerType$n[16]+tabPerType$n[20]+tabPerType$n[23]))
  tmpNorm[j,2] <- c
}
tmpNorm <- as.data.frame(tmpNorm)
names(tmpNorm)[1] <- paste("cluster")
names(tmpNorm)[2] <- paste("freq")
lungnormalClus <- tmpNorm

lungcoreClus <- rbind(lung1coreClus, lung2coreClus, lung3coreClus, lung4coreClus, lung5coreClus, lung6coreClus)
lungcoreClus <- select(lungcoreClus, one_of(c("cluster", "n")))
tmpCore <- matrix(ncol=2, nrow=dim(tabPerClus)[1])
for(j in 1:dim(tabPerClus)[1]){
  iter <- as.numeric(j)
  tmpCore[j,1] <- iter
  a <- lapply(filter(lungcoreClus, lungcoreClus$cluster == j) %>% summarise_at(c("n"),sum,na.rm=TRUE), as.numeric)
  b <- a[[1]]
  c <- (b/(tabPerType$n[1]+tabPerType$n[5]+tabPerType$n[9]+tabPerType$n[13]+tabPerType$n[17]+tabPerType$n[21]))
  tmpCore[j,2] <- c
}
tmpCore <- as.data.frame(tmpCore)
names(tmpCore)[1] <- paste("cluster")
names(tmpCore)[2] <- paste("freq")
lungcoreClus <- tmpCore

lungmiddleClus <- rbind(lung1middleClus,lung2middleClus,lung3middleClus,lung4middleClus,lung5middleClus)
lungmiddleClus <- select(lungmiddleClus, one_of(c("cluster", "n")))
tmpMiddle <- matrix(ncol=2, nrow=dim(tabPerClus)[1])
for(j in 1:dim(tabPerClus)[1]){
  iter <- as.numeric(j)
  tmpMiddle[j,1] <- iter
  a <- lapply(filter(lungmiddleClus, lungmiddleClus$cluster==j) %>% summarise_at(c("n"),sum,na.rm=TRUE), as.numeric)
  b <- a[[1]]
  c <- (b/(tabPerType$n[3]+tabPerType$n[7]+tabPerType$n[11]+tabPerType$n[15]+tabPerType$n[19]))
  tmpMiddle[j,2] <- c
}
tmpMiddle <- as.data.frame(tmpMiddle)
names(tmpMiddle)[1] <- paste("cluster")
names(tmpMiddle)[2] <- paste("freq")
lungmiddleClus <- tmpMiddle

lungedgeClus <- rbind(lung1edgeClus,lung2edgeClus,lung3edgeClus,lung4edgeClus,lung5edgeClus,lung6edgeClus)
lungedgeClus <- select(lungedgeClus, one_of(c("cluster", "n")))
tmpEdge <- matrix(ncol=2, nrow=dim(tabPerClus)[1])
for(j in 1:dim(tabPerClus)[1]){
  iter <- as.numeric(j)
  tmpEdge[j,1] <- iter
  a <- lapply(filter(lungedgeClus, lungedgeClus$cluster == j) %>% summarise_at(c("n"),sum,na.rm=TRUE), as.numeric)
  b <- a[[1]]
  c <- (b/(tabPerType$n[2]+tabPerType$n[6]+tabPerType$n[10]+tabPerType$n[14]+tabPerType$n[18]+tabPerType$n[22]))
  tmpEdge[j,2] <- c
}
tmpEdge <- as.data.frame(tmpEdge)
names(tmpEdge)[1] <- paste("cluster")
names(tmpEdge)[2] <- paste("freq")
lungedgeClus <- tmpEdge

# function to calculate diversity score for a given q
calcqD <- function(dataset, q){
  diversity <- 0.0;
  for(row in 1:dim(dataset)[1]){
    diversity <- diversity + (dataset$freq[row])^q
  }
  diversity <- diversity^(1/(1-q))
}

# calculates diversity scores based on cluster frequencies for a range of q
qRange <- logspace(-2, 2, n = 1000)
qDlung1normal <- calcqD(lung1normalClus,qRange)
qDlung1normaltab <- as.data.frame(cbind(qRange,qDlung1normal))
qDlung1core <- calcqD(lung1coreClus,qRange)
qDlung1coretab <- as.data.frame(cbind(qRange,qDlung1core))
qDlung1middle <- calcqD(lung1middleClus,qRange)
qDlung1middletab <- as.data.frame(cbind(qRange,qDlung1middle))
qDlung1edge <- calcqD(lung1edgeClus,qRange)
qDlung1edgetab <- as.data.frame(cbind(qRange,qDlung1edge))

qDlung2normal <- calcqD(lung2normalClus,qRange)
qDlung2normaltab <- as.data.frame(cbind(qRange,qDlung2normal))
qDlung2core <- calcqD(lung2coreClus,qRange)
qDlung2coretab <- as.data.frame(cbind(qRange,qDlung2core))
qDlung2middle <- calcqD(lung2middleClus,qRange)
qDlung2middletab <- as.data.frame(cbind(qRange,qDlung2middle))
qDlung2edge <- calcqD(lung2edgeClus,qRange)
qDlung2edgetab <- as.data.frame(cbind(qRange,qDlung2edge))

qDlung3normal <- calcqD(lung3normalClus,qRange)
qDlung3normaltab <- as.data.frame(cbind(qRange,qDlung3normal))
qDlung3core <- calcqD(lung3coreClus,qRange)
qDlung3coretab <- as.data.frame(cbind(qRange,qDlung3core))
qDlung3middle <- calcqD(lung3middleClus,qRange)
qDlung3middletab <- as.data.frame(cbind(qRange,qDlung3middle))
qDlung3edge <- calcqD(lung3edgeClus,qRange)
qDlung3edgetab <- as.data.frame(cbind(qRange,qDlung3edge))

qDlung4normal <- calcqD(lung4normalClus,qRange)
qDlung4normaltab <- as.data.frame(cbind(qRange,qDlung4normal))
qDlung4core <- calcqD(lung4coreClus,qRange)
qDlung4coretab <- as.data.frame(cbind(qRange,qDlung4core))
qDlung4middle <- calcqD(lung4middleClus,qRange)
qDlung4middletab <- as.data.frame(cbind(qRange,qDlung4middle))
qDlung4edge <- calcqD(lung4edgeClus,qRange)
qDlung4edgetab <- as.data.frame(cbind(qRange,qDlung4edge))

qDlung5normal <- calcqD(lung5normalClus,qRange)
qDlung5normaltab <- as.data.frame(cbind(qRange,qDlung5normal))
qDlung5core <- calcqD(lung5coreClus,qRange)
qDlung5coretab <- as.data.frame(cbind(qRange,qDlung5core))
qDlung5middle <- calcqD(lung5middleClus,qRange)
qDlung5middletab <- as.data.frame(cbind(qRange,qDlung5middle))
qDlung5edge <- calcqD(lung5edgeClus,qRange)
qDlung5edgetab <- as.data.frame(cbind(qRange,qDlung5edge))

qDlung6normal <- calcqD(lung6normalClus,qRange)
qDlung6normaltab <- as.data.frame(cbind(qRange,qDlung6normal))
qDlung6core <- calcqD(lung6coreClus,qRange)
qDlung6coretab <- as.data.frame(cbind(qRange,qDlung6core))
qDlung6edge <- calcqD(lung6edgeClus,qRange)
qDlung6edgetab <- as.data.frame(cbind(qRange,qDlung6edge))

# join datasets together for plotting aggregated dataset
qDlungnormal <- calcqD(lungnormalClus,qRange)
qDlungnormaltab <- as.data.frame(cbind(qRange,qDlungnormal))
qDlungcore <- calcqD(lungcoreClus,qRange)
qDlungcoretab <- as.data.frame(cbind(qRange,qDlungcore))
qDlungmiddle <- calcqD(lungmiddleClus,qRange)
qDlungmiddletab <- as.data.frame(cbind(qRange,qDlungmiddle))
qDlungedge <- calcqD(lungedgeClus,qRange)
qDlungedgetab <- as.data.frame(cbind(qRange,qDlungedge))

# plots all data together
ggplot(data=qDlung1normaltab, aes(x=qRange, y=qDlung1normal, color="normal")) +
  geom_smooth(size=6) + 
  scale_x_log10() +
  geom_smooth(data=qDlung1coretab, aes(x=qRange, y=qDlung1core, color="core"),size=6) +
  geom_smooth(data=qDlung1middletab, aes(x=qRange, y=qDlung1middle, color="middle"),size=6) +
  geom_smooth(data=qDlung1edgetab, aes(x=qRange, y=qDlung1edge, color="edge"),size=6) +
  geom_smooth(data=qDlung2normaltab, aes(x=qRange, y=qDlung2normal, color="normal"),size=6) +
  geom_smooth(data=qDlung2coretab, aes(x=qRange, y=qDlung2core, color="core"),size=6) +
  geom_smooth(data=qDlung2middletab, aes(x=qRange, y=qDlung2middle, color="middle"),size=6) +
  geom_smooth(data=qDlung2edgetab, aes(x=qRange, y=qDlung2edge, color="edge"),size=6) +
  geom_smooth(data=qDlung3normaltab, aes(x=qRange, y=qDlung3normal, color="normal"),size=6) +
  geom_smooth(data=qDlung3coretab, aes(x=qRange, y=qDlung3core, color="core"),size=6) +
  geom_smooth(data=qDlung3middletab, aes(x=qRange, y=qDlung3middle, color="middle"),size=6) +
  geom_smooth(data=qDlung3edgetab, aes(x=qRange, y=qDlung3edge, color="edge"),size=6) +
  geom_smooth(data=qDlung4normaltab, aes(x=qRange, y=qDlung4normal, color="normal"),size=6) +
  geom_smooth(data=qDlung4coretab, aes(x=qRange, y=qDlung4core, color="core"),size=6) +
  geom_smooth(data=qDlung4middletab, aes(x=qRange, y=qDlung4middle, color="middle"),size=6) +
  geom_smooth(data=qDlung4edgetab, aes(x=qRange, y=qDlung4edge, color="edge"),size=6) +
  geom_smooth(data=qDlung5normaltab, aes(x=qRange, y=qDlung5normal, color="normal"),size=6) +
  geom_smooth(data=qDlung5coretab, aes(x=qRange, y=qDlung5core, color="core"),size=6) +
  geom_smooth(data=qDlung5middletab, aes(x=qRange, y=qDlung5middle, color="middle"),size=6) +
  geom_smooth(data=qDlung5edgetab, aes(x=qRange, y=qDlung5edge, color="edge"),size=6) +
  geom_smooth(data=qDlung6normaltab, aes(x=qRange, y=qDlung6normal, color="normal"),size=6) +
  geom_smooth(data=qDlung6coretab, aes(x=qRange, y=qDlung6core, color="core"),size=6) +
  geom_smooth(data=qDlung6edgetab, aes(x=qRange, y=qDlung6edge, color="edge"),size=6) +
  labs(title="Diversity Score",x="q", y = "qD") +
  scale_color_manual(name="", 
                     breaks = c("core","edge","middle","normal"),
                     values = c("darkgreen","green3","forestgreen","black"),
                     guide = guide_legend(override.aes=aes(fill=NA))) +
  theme(rect = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))

# Plotting joined dataset
ggplot(data=qDlungnormaltab, aes(x=qRange, y=qDlungnormal)) +
   geom_smooth(size=6, color="black") + #color="dark gray",
   scale_x_log10() +
   geom_smooth(data=qDlungcoretab, aes(x=qRange, y=qDlungcore),size=6, color="darkgreen") +
   geom_smooth(data=qDlungmiddletab, aes(x=qRange, y=qDlungmiddle),size=6, color="forestgreen") +
   geom_smooth(data=qDlungedgetab, aes(x=qRange, y=qDlungedge),size=6, color="green3") +
   labs(title="Diversity Score",x="q", y = "qD") +
   ylim(0,25) +
   theme(rect = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"))
