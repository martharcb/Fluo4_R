
#Creating subsets of data for graphing individual traces from calcium imaging data
#Last updated on 1-6-20 by Erik Larsen and Martha Bhattacharya

#Step 1: Load libraries
library(tidyverse)
library(ggplot2)
library(purrr)
library(reshape2)
library(dplyr)


# Import the fluorescence value CSVs.
# Make sure the csv's are trimmed appropriately. Delete extra rows so that the first row containing text is the first row of the file.
B1 = read.csv("N:/Martha/Ca_imaging/2018.12.14.il31.cap.wt.mut/wt1.tifs/stats.csv", header =F)

# Create frame (t) and time (sec) vectors

#     you should set your time vecLB1 <- nrow(B1)
t <- 1:LB1
tsec <- t*2
FLUOwt1 = B1
#normalize the data
for (i in 1:ncol(FLUOwt1)){
  FLUOwt1[,i] = (FLUOwt1[,i]) / (FLUOwt1[2:7,i]) 
}

#RERUN FOR WT2
B1 = read.csv("N:/Martha/Ca_imaging/2018.12.14.il31.cap.wt.mut/wt2 tifs/stats.csv", header =F)

# Create frame (t) and time (sec) vectors

#     you should set your time vecLB1 <- nrow(B1)
t <- 1:LB1
tsec <- t*2
FLUOwt2 = B1
#normalize the data
for (i in 1:ncol(FLUOwt2)){
  FLUOwt2[,i] = (FLUOwt2[,i]) / (FLUOwt2[2:7,i]) 
}

#run the code above for each data file you need to import
#You will need to crop files with longer hiK treatments by using something like the following after import:
     FLUOwt2 = FLUOwt2[1:127,]
  tsec = tsec[1:127]
  
  #RERUN FOR MUT1
  B1 = read.csv("N:/Martha/Ca_imaging/2018.12.14.il31.cap.wt.mut/mu1 tifs/stats.csv", header =F)
  
  # Create frame (t) and time (sec) vectors
  
  #     you should set your time vecLB1 <- nrow(B1)
  t <- 1:LB1
  tsec <- t*3
  FLUOmu1 = B1
  #normalize the data
  for (i in 1:ncol(FLUOmu1)){
    FLUOmu1[,i] = (FLUOmu1[,i]) / (FLUOmu1[2:7,i]) 
  }
  
  #run the code above for each data file you need to import
  #You will need to crop files with longer hiK treatments by using something like the following after import:
  FLUOmu1 = FLUOmu1[1:127,]
  tsec = tsec[1:127]
  
  #RERUN FOR MUT2
  B1 = read.csv("N:/Martha/Ca_imaging/2018.12.14.il31.cap.wt.mut/mu2 tifs/stats.csv", header =F)
  
  # Create frame (t) and time (sec) vectors
  
  #     you should set your time vecLB1 <- nrow(B1)
  t <- 1:LB1
  tsec <- t*3
  FLUOmu2 = B1
  #normalize the data
  for (i in 1:ncol(FLUOmu2)){
    FLUOmu2[,i] = (FLUOmu2[,i]) / (FLUOmu2[2:7,i]) 
  }
  
  #run the code above for each data file you need to import
  #You will need to crop files with longer hiK treatments by using something like the following after import:
  FLUOmu2 = FLUOmu2[1:127,]
  tsec = tsec[1:127]
  
  
  
#first create a data frame with the columns containing the cells you want to graph and the time vector


FLUOsubset = as.data.frame(cbind(tsec, FLUOmu1$V311, FLUOmu1$V403, FLUOmu1$V446, FLUOmu2$V2, FLUOmu2$V10, FLUOmu2$V21, FLUOmu2$V61, FLUOmu2$V81, FLUOmu2$V83, FLUOmu2$V110, FLUOmu2$V134, FLUOmu2$V141, FLUOmu2$V143, FLUOmu2$V146, FLUOmu2$V147, FLUOmu2$V152, FLUOmu2$V178, FLUOmu2$V225, FLUOmu2$V270, FLUOmu2$V302, FLUOmu2$V322, FLUOmu2$V326, FLUOmu2$V389, FLUOmu2$V402))

FLUO_WTsubset = as.data.frame(cbind(tsec, FLUOwt1$V63, FLUOwt1$V175, FLUOwt1$V241, FLUOwt2$V7, FLUOwt2$V54, FLUOwt2$V60, FLUOwt2$V92, FLUOwt2$V111, FLUOwt2$V122, FLUOwt2$V162, FLUOwt2$V214, FLUOwt2$V230, FLUOwt2$V236, FLUOwt2$V274, FLUOwt2$V275, FLUOwt2$V284, FLUOwt2$V288, FLUOwt2$V296, FLUOwt2$V302, FLUOwt2$V323, FLUOwt2$V326))

#colnames(FLUOsubset) = c("tsec", "mu1227", "mu1248", "mu1272", "wt163", "wt2302", "wt2326")
mutmean = rowMeans(FLUOsubset[1:nrow(FLUOsubset),2:ncol(FLUOsubset)])
wtmean = rowMeans(FLUO_WTsubset[1:nrow(FLUO_WTsubset),2:ncol(FLUO_WTsubset)])
mutgraph = cbind(tsec, mutmean, wtmean)

mutgraph = as.data.frame(mutgraph)
CaImaging = ggplot(data = mutgraph[2:127,], aes(x = mutgraph[2:127,1])) +
  theme_classic() +
  #trace 1
  geom_line(aes(y=mutgraph[2:127,2], color="red")) +
  geom_point(aes(y=mutgraph[2:127,2], color="red")) +
  #trace 2
  geom_line(aes(y=mutgraph[2:127,3], color="orange")) +
  geom_point(aes(y=mutgraph[2:127,3], color="orange")) +
  # #trace 3
  # geom_line(aes(y=FLUOsubset$mu1272[2:127], color="yellow")) + 
  # geom_point(aes(y=FLUOsubset$mu1272[2:127], color="yellow")) +
  # #trace 4
  # geom_line(aes(y=FLUOsubset$wt163[2:127], color="blue")) +
  # geom_point(aes(y=FLUOsubset$wt163[2:127], color="blue"))+
  # #trace 5
  # geom_line(aes(y=FLUOsubset$wt2302[2:127], color="green")) +
  # geom_point(aes(y=FLUOsubset$wt2302[2:127], color="green")) +
  # #trace 6
  # geom_line(aes(y=FLUOsubset$wt2326[2:127], color="purple")) + 
  # geom_point(aes(y=FLUOsubset$wt2326[2:127], color="purple")) +
  #add guidelines to the graph, if desired:
  geom_vline(xintercept = 27, color = "forest green") +
  geom_vline(xintercept = 270, color = "red") +
  geom_vline(xintercept = 330, color = "purple") + 
  #+ geom_hline(yintercept = 1.0, color = "black", linetype = "solid")
  coord_cartesian(xlim = c(3,180), ylim = c(0.5,3)) +
  labs(x = "Time (s)", y = "Normalized Fluorescence (a.u.)") +
  theme(plot.title = element_text(hjust = 0.5))

CaImaging
#export individual cell plots
mutgraph = as.data.frame(cbind(tsec, wtmean, mutmean, FLUOmu1$V446, FLUOmu2$V2, FLUOmu2$V81))
mutgraph = as.data.frame(cbind(mutgraph, FLUOwt2$V274, FLUOwt2$V302, FLUOwt2$V323, FLUOwt2$V54))
#FLUOwt1$V241, FLUOwt1$V274, FLUOwt1$V302, FLUOwt1$V323, FLUOwt1$V54))
write_csv(as.data.frame(mutgraph), path = "N:\\Martha\\Ca_imaging\\2018.12.14.il31.cap.wt.mut\\meangraph7.csv")

#################################################################################################################
#get avg without Wt2-275 and 288, which have jumps in fluorescence that change WT 
FLUO_WTsubset2 = as.data.frame(cbind(tsec, FLUOwt1$V63, FLUOwt1$V175, FLUOwt1$V241, FLUOwt2$V7, FLUOwt2$V54, FLUOwt2$V60, FLUOwt2$V92, FLUOwt2$V111, FLUOwt2$V122, FLUOwt2$V162, FLUOwt2$V214, FLUOwt2$V230, FLUOwt2$V236, FLUOwt2$V274, FLUOwt2$V284, FLUOwt2$V296, FLUOwt2$V302, FLUOwt2$V323, FLUOwt2$V326))

wtmean2 = rowMeans(FLUO_WTsubset2[1:nrow(FLUO_WTsubset2),2:ncol(FLUO_WTsubset2)])
mutgraph2 = cbind(tsec, mutmean, wtmean2)
mutgraph2 = as.data.frame(mutgraph2)

CaImaging = ggplot(data = mutgraph2[2:127,], aes(x = mutgraph2[2:127,1])) +
  theme_classic() +
  #trace 1
  geom_line(aes(y=mutgraph2[2:127,2], color="red")) +
  geom_point(aes(y=mutgraph2[2:127,2], color="red")) +
  #trace 2
  geom_line(aes(y=mutgraph2[2:127,3], color="orange")) +
  geom_point(aes(y=mutgraph2[2:127,3], color="orange")) +
#add guidelines to the graph, if desired:
geom_vline(xintercept = 27, color = "forest green") +
  geom_vline(xintercept = 270, color = "red") +
  geom_vline(xintercept = 330, color = "purple") + 
  #+ geom_hline(yintercept = 1.0, color = "black", linetype = "solid")
  coord_cartesian(xlim = c(3,180), ylim = c(0.5,3)) +
  labs(x = "Time (s)", y = "Normalized Fluorescence (a.u.)") +
  theme(plot.title = element_text(hjust = 0.5))

write_csv(as.data.frame(mutgraph2), path = "N:\\PAPER ASSEMBLY\\Itch paper\\Calcium Imaging Plots\\meangraph8.csv")
write_csv(as.data.frame(FLUOwt2$V214), path = "N:\\PAPER ASSEMBLY\\Itch paper\\Calcium Imaging Plots\\meangraph10.csv")

##plotting individual cells with ggplot
FLUOwt1 = as.data.frame(cbind(tsec, FLUOwt1))
cellplot = ggplot(data = FLUOwt1, aes(x = FLUOwt1[,1])) + geom_point(aes(y=FLUOwt1[,115], color = "red")) + geom_line(aes(y=FLUOwt1[,115], color = "red")) + coord_cartesian(xlim = c(3,180), ylim = c(0.5,10)) + labs(x = "Time (s)", y = "Normalized Fluorescence (a.u.)") + theme(plot.title = element_text(hjust = 0.5))
