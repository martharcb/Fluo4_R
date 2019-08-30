#This program will take logs of cell fluorescence values over time from a calcium
#imaging experiment and use it to determine the number of responders and their response magnitudes
#Updated by Martha Bhattacharya on 8/30/19

#set working directory to server folder
setwd("M:/Calcium Imaging Analysis")

#load packages needed for the program to run. Current locations are Martha's 
#hard drive but should be changed by the user.
library("tidyverse", lib.loc="M:/Calcium Imaging Analysis/packagesR")
library("tcltk", lib.loc="M:/Calcium Imaging Analysis/packagesR")
library("ggplot2", lib.loc="M:/Calcium Imaging Analysis/packagesR")

#make variables for data directory and output directory
data.dir <- "M:/Calcium Imaging Analysis/rawdata"
results.dir <- "M:/Calcium Imaging Analysis/outputR"

#import and make data frame and time vectors
well7 <- read.csv("well7_fluo4_Results.csv")
well7df <- data.frame(well7)
t <- well7df[,1]
tsec <- (t*3)-3

#import bkgd files and subtract background
well7bkgd3 <- read.csv("well7bkgd3.csv")
avg_bk <- rowMeans(well7bkgd3[,2:4])
fluo <- well7df[,2:ncol(well7df)]-avg_bk

#determine baseline and normalize each cell to its baseline
ratio <- data.frame(fluo)
ratio <- sweep(fluo,MARGIN=2,FUN = "/",STATS=colMeans(fluo[2:6,]))

#put all required elements into a single data table
ratio_t <- data.frame(cbind(t,tsec,ratio))

#extract information from each cell
#first create a list of cell numbers
cellnum <- as.numeric(gsub("Mean", "", colnames(ratio)))
#then get peaks from time windows
peak1 <- vector("double", ncol(ratio))
cap1 <- vector("double",ncol(ratio))
hik1 <- vector("double",ncol(ratio))
for (i in seq_along(ratio)) {
    peak1[i] <- max(ratio[11:30,i])
    cap1[i] <- max(ratio[71:85,i])
    hik1[i] <- max(ratio[91:105,i])
}
#put all dta together in a single frame
allpeaks <- data.frame(cbind(cellnum,peak1,cap1,hik1))


#making individual plots as separate png files 

dir.create("plots5")
colNames <- names(ratio_t)
for(i in colNames){
  outpath = file.path(getwd(), "plots5", paste("cell_",i,".png", sep=""))
  ggplot(ratio_t, aes_string(x=t, y = i), width=4, height=3)+geom_point()+geom_line()+coord_cartesian(ylim=c(0.8,2))
  ggsave(filename = outpath,plot = last_plot(), width=4, height=3)
}
