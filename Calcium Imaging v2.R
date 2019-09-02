
# set the current directory with setwd()
setwd("M:/Erik/Data/Calcium Imaging Analysis/Batch 2")

install.packages("tidyverse")
install.packages("reshape2")
install.packages("stats")
install.packages("caret")

library(ggplot2)
library(plyr)
library(ggrepel)
library(calibration)
library(reshape2)
library(gridExtra)
library(stats)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(purrr)
library(broom)
library(stringr)



##### Getting the data #####
# Make a vector (column) of all the csv files in the newly set working directory
csvs = dir(pattern = "\\.csv$")


# Import the csv files, each as their own dataframes

# File 1 is GREEN; file 3 is GREEN BKGD
# File 2 is RED; file 4 is RED BKGD
green = read.csv(csvs[1])
greenbkgd = read.csv(csvs[3])

# Make a new dataframe with normalized signal

# *****Analyze file by file.***** In this data set, the first frame seems to be an artefact, so will be skipped.
green = green[2:130,3:ncol(green)]
greenbkgd = greenbkgd[2:130,3]


gsub = green[,] - greenbkgd

# Apply for red
red = read.csv(csvs[2])

red = red[2:130,3:ncol(red)]
redbkgd = read.csv(csvs[4])
redbkgd = redbkgd[2:130,3]

rsub = red[,] - redbkgd


# Make a time vector (s) from the original GREEN dataframe (should not matter: GREEN and RED should have the same number of rows)
t = 1:nrow(green)

# Create a new dataframe, containing the ratioio of activity normalized to background
ratio = gsub/rsub


# Plot the first cell over time as a check
plot(t, ratio[,1], xlim = c(0,170), ylim = c(-1,4))



##### Prepping the data for creating plots of cellular responses #####

# Copy the response dataframe
Adjust = ratio


# Adjust the baseline fluorescence for each cell to a value of 1 at time 0 in the new response dataframe
for (i in 1:ncol(ratio)){
  for (j in 1:nrow(ratio)){
    Adjust[j,i] = if_else(ratio[1,i] < 1, ratio[j,i]+(1-ratio[1,i]), ratio[j,i]+(1-ratio[1,i])) 
  }
}



# Find the minimum, maximum for plotting scales and removing artefacts
min(Adjust[,])
max(Adjust[,])
which(Adjust == min(Adjust), arr.ind = TRUE)
which(Adjust == max(Adjust), arr.ind = TRUE)

# Manual attempt at interpolation
for (i in 1:length(Adjust)){
  for (j in 1:nrow(Adjust)) {
      Adjust[j,i] = if_else(Adjust[j,i] > 5, ((Adjust[j,i+1] + Adjust[j,i-1])/2),
                            if_else(Adjust[j,i] < 0, ((Adjust[j,i+1] + Adjust[j,i-1])/2), Adjust[j,i]))
  }
} 

Adjust[Adjust > 5 | Adjust < 0] = Adjust[Adjust]


Adjust %>%
      mutate_all(Adjust())


                      for (i = 1) {
                        ((Adjust[j,i] + Adjust[j,i+1])/2)
    } 
                      for (i = 175) {
                        (Adjust[j,i] + Adjust[j,i-1])/2
    }
      Adjust[j,i] = (Adjust[j,i+1] + Adjust[j,i-1])/2, Adjust[j,i] = Adjust[j,i])
  }
}






##### First Agonist Threshold for Classification #####




# Set an arbitrary baseline average FL (pre-agonist) for all cells prior to addition of the first agonist in the 10th frame
# Take the average fluorescence over the first few frames
ratio_avg_baseline = (colSums(Adjust[1:9, 1:length(Adjust)], dims = 1))/9


# Calculate the response threshold for the first agonist (at frame 25; 15 frames after exposure)
First_agonist_threshold = ratio_avg_baseline*0.10 + ratio_avg_baseline


# Determine responders to the first agonist
# Create an empty vector to fill
First_agonist_responses = vector()

# Fill the vector with cells that responded to the agonist or not, based on the pre-determined, arbitrary threshold of what constitutes a "response"
for (i in 1:length(First_agonist_threshold)){
  First_agonist_responses[i] = if_else(Adjust[25,i] > First_agonist_threshold[i], "Responder", "Non-responder")
}


Adjust[11,38] - ratio_avg_baseline[38]


# See how many and also which cells are classified as "responders"
length(First_agonist_responses[First_agonist_responses == "Responder"])
which(First_agonist_responses == "Responder")






##### Capsaiscin Threshold for Classification #####



# Determine responders to capsaicin
# Establish a baseline, accounting for the mean fluorescence in the 8 frames leading up to the Capsaicin addition
Capsaicin_baseline = (colSums(Adjust[80:87,1:length(Adjust)-1], dims = 1))/8


# Set the threshold at frame 91; 2 frames after exposure
Capsaicin_threshold = Capsaicin_baseline*0.05 + Adjust[91,1:length(Adjust)-1]


# Create an empty vector to fill
Capsaicin_responses = vector()


# Make response calls after addition in the 89th frame
for (i in 1:length(Capsaicin_threshold)){
  Capsaicin_responses[i] = if_else(Adjust[91,i] - Capsaicin_baseline[i] > 0.05, "Responder", "Non-responder")
}

# See how many and also which cells are classified as "responders"
length(Capsaicin_responses[Capsaicin_responses == "Responder"])
which(Capsaicin_responses == "Responder")




##### Neural Threshold for Classification #####


# Determine neurons
# Establish a baseline, accounting for the mean in the 8 frames leading up to the High K addition
Neural_baseline = (colSums(Adjust[100:107,1:length(Adjust)-1], dims = 1))/8


# Set the threshold at 5% above the baseline established 8 frames before High K addition
Neural_threshold = Neural_baseline*0.05 + Adjust[111,1:length(Adjust)-1]


# Create an empty vector to fill
High_K_responses = vector()

# Make response calls in the 110th frame to High K addition in the 109th frame
for (i in 1:length(Neural_threshold)){
  High_K_responses[i] = if_else(Capsaicin_responses[i] == "Responder", "Neuron", 
                                if_else(Adjust[110,i] - Neural_baseline[i] > 0.5, "Neuron", "Non-neuron"))
}

High_K_duration = vector()

for (i in 1:length(Neural_threshold)){
  for (j in 1:length(t)){
    for (k in 110:j) {
      if ((Adjust[(110:j),i]) > Adjust[110,i]) {
        High_K_duration[i] = length(k)}
    }
  }
}

    High_K_duration[i] = 
      if_else(Adjust[(110:j),i]) > Adjust[110,i], length(t), 
  }
  
}

# See how many and also which cells are classified as "Neurons"
length(High_K_responses[High_K_responses == "Neuron"])
which(High_K_responses == "Neuron")


dummy = vector()
for (i in 1:length(High_K_responses)){
  if (High_K_responses[i] == "Neuron") {
    dummy[i] = Adjust[11,i]
  } else {
    dummy[i] = 0
  }
}

length(dummy[dummy] != 0)
which(dummy == 0)

which(High_K_responses == "Non-neuron" & First_agonist_responses == "Responder")


##### Create a new dataframe that will store binary character response classifications of each cell to each reagent along with the response magnitude and duration #####


# Create a dataframe containing the responses of each cell for each agonist





Responses = data.frame()

for (i in 1:length(Adjust)) {
  # 1st Ag Mag
  Responses[i,1] = Adjust[11,i]
  # 2nd Ag Mag
  Responses[i,2] = Adjust [91,i]
  # High K Mag
  Responses[i,3] = Adjust [111,i]
  # 1st Ag Duration
  Responses[i,4] = Adjust[(c(11:13)/3),1]i]
  # 2nd Ag Duration
  Responses[i,5] = Adjust[,i]
  # High K Duration
  Responses[i,6] = Adjust[,i]
}





Response_DF = data.frame(First_agonist = c(First_agonist_responses), First_ag_Mag = Responses[,1], Capsaicin = c(Capsaicin_responses), Capsaicin_Mag = Responses[,2], High_K = c(High_K_responses), High_K_Mag = Responses[,3])





x = Response_DF[,c(2,4,6)]
y = Response_DF[,c(1,3,5)]

par(mfrow = c(1,3))
boxplot(x[,], main = names(Response_DF[2,4,6]))



par(mfrow = c(1,3))
for (i in 1:3) {
  boxplot(x[,i], main =names(Response_DF[,c(2,4,6)])[i], ylim = c(-2,3))
}



featurePlot(x=x, y=y, plot = "ellipse")
  




# For calculations
dummy = Adjust[10,11] - ratio_avg_baseline[11]


# For outliers
artefacts = Adjust

# Revise outliers
for (i in 1:ncol(Adjust[,ncol(Adjust)-1])){
  for (j in 1:nrow(Adjust)){
    artefacts[j,i] = if_else(Adjust[j,i] > 2, Adjust[j,i]+(1-Adjust[1,i]), Adjust[j,i]+(1-Adjust[1,i])) 
  }
}

Adjust = cbind(t, Adjust)

# GGplot of a cell
CaImaging = ggplot(data = Adjust, aes(x = Adjust$t, y = Adjust[,39])) + geom_point(alpha = 0.3, color = "red") + coord_cartesian(xlim = c(0,130), ylim = c(0.5,2.5)) + labs(x = "Time (Frame#; 4s/Frame)", y = "Normalized Fluorescence (a.u.)", title = "In vitro WT DRG Calcium Imaging Responses") + theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept = 9, color = "forest green", linetype = "dashed") + geom_vline(xintercept = 89, color = "orange", linetype = "dashed") + geom_vline(xintercept = 109, color = "navy", linetype = "dashed") + geom_hline(yintercept = 1, color = "black", linetype = "solid") + theme_light()

CaImaging





##### Save quick plots of each cell into a subdirectory on the server

# Graph each cell, put them in a subdirectory called "plots"
# graph <- function(ratio) {         this isn't working yet!
  dir.create("plots")
  for (i in seq_along(Adjust[,ncol(Adjust)-1])) {
    #tempname = paste("cell_",i,".png", sep="")
    outpath = file.path("M:\\Erik\\Data\\Calcium Imaging Analysis\\Batch 2\\plots", paste("cell_",i,".png", sep=""))
    png(filename = outpath) 
    plot(t, Adjust[,i], xlim = c(0,145), ylim = c(-5,1300))
    #ggsave(filename=outpath)
    dev.off()
  }
#}

##### Save GGplots of each cell into a subdirectory on the server #####
  ## For looping through ggplots ##
  for (i in seq_along(Adjust[,ncol(Adjust)-1])) {
    CaImaging[i] = ggplot(data = Adjust, aes(x = Adjust$t, y = Adjust[,i])) + geom_point(alpha = 0.2, color = "red") + coord_cartesian(xlim = c(0,145), ylim = c(0,3)) + labs(x = "Time (s*4)", y = "Normalized Fluorescence", title = "In-vitro DRG Calcium Imaging Responses to Various Agonists") + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0))
  }
  
ROIs = unique(Adjust)
  
  
  dir.create("plots2")
  for (i in seq_along(Adjust[,ncol(Adjust)-1])) {
    #tempname = paste("cell_",i,".png", sep="")
    outpath = file.path("M:\\Erik\\Data\\Calcium Imaging Analysis\\Batch 2\\plots2", paste("cell_",i,".png", sep=""))
    CaImaging = ggplot(data = Adjust, aes(x = Adjust$t, y = Adjust[,i])) + geom_point(alpha = 0.2, color = "red") + coord_cartesian(xlim = c(0,145), ylim = c(0,3)) + labs(x = "Time (s*4)", y = "Normalized Fluorescence of Cell_", paste(i), title = "In-vitro DRG Calcium Imaging Responses to Various Agonists") + theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0))
    print(CaImaging)
    ggsave(outpath[i], CaImaging[i], device = "png", scale = 0.5)
    #ggsave(filename=outpath)
    dev.off()
  }
  
  

  

##### K Means Clustering #####  

  
  
kclust = kmeans(Response_DF, centers = 3, nstart = 3)

summary(kclust)

kclusts = tibble(k = 1:9) %>%
  mutate(
    kclust = map(k, ~kmeans(Response_DF, x = 9)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, ratio[,ncol(Adjust)-1])
  )

clusters <- kclusts %>%
  unnest(tidied)

assignments <- kclusts %>% 
  unnest(augmented)

clusterings <- kclusts %>%
  unnest(glanced, .drop = TRUE)


p1 <- ggplot(assignments, aes(, )) +
  geom_point(aes(color = .cluster)) + 
  facet_wrap(~ k)
p1








##### Notes/Garbage code #####
# Create the dataframe, with the three columns to be used to classify neurons (subpopulations) by their responses to given agonists and the magnitude of those responses (Agonist, (binary response/no response; neuron/non-neuron), Magnitude of response)
Response_DF = data.frame(Agonist=character(length(Adjust)-1), Response=character(3*length(Adjust)-3), Magnitude= numeric(3*length(Adjust)-3), stringsAsFactors = FALSE)


# Insert agonist type into the dataframe
# 1st 175 = 1st agonist
Response_DF$Agonist[1:(length(Adjust)-1)] = "First Agonist"
# 2nd 175 = capsaicin
Response_DF$Agonist[length(Adjust):(2*length(Adjust)-1)] = "Capsaicin"
# 3rd 175 = high potassium
Response_DF$Agonist[(2*length(Adjust)):nrow(Response_DF)] = "High K"


# Add in whether they responded according to pre-determined thresholds
# 1st agonist responses
Response_DF$Response[1:(length(Adjust)-1)] = c(First_agonist_responses)
# capsaicin responses
Response_DF$Response[length(Adjust):(2*length(Adjust)-1)] = c(Capsaicin_responses)
# high potassium responses
Response_DF$Response[(2*length(Adjust)):nrow(Response_DF)] = c(High_K_responses)


# Insert the magnitudes of fluorescence of the responders ("responders", "responders", "neurons"). Otherwise, 0.

for (i in 1:length(Response_DF$Agonist)){
  if (Response_DF$Agonist[i] == "First Agonist")
    if (Response_DF$Response[i] == "Responder") {
      Response_DF$Magnitude[i] = Adjust[11,i] 
    } else if (Response_DF$Agonist[i] == "Capsaicin")
      if (Response_DF$Response[i] == "Responder"){
        Response_DF$Magnitude[i] = Adjust[91,i]
      } else if (Response_DF$Agonist[i] == "High K")
        if (Response_DF$Response[i] == "Neuron") {
          Response_DF$Magnitude[i] = Adjust[110,i]
        } else if (Response_DF[i] == "First Agonist")
          if (Response_DF$Response[i] != "Responder") {
            Response_DF$Magnitude[i] = 0
          } else if (Response_DF$Agonist[i] == "Capsaicin")
            if (Response_DF$Response[i] != "Responder") {
              Response_DF$Magnitude[i] = 0
            } else {
              Response_DF$Magnitude[i] = 0
            }
}

Response_DF$Magnitude[i] = Adjust[11,i] {
  else if (Response_DF$Response[i] == "")
}
for (j in length(Adjust):(2*length(Adjust)-1)){
  Response_DF$Magnitude[length(Adjust):(2*length(Adjust)-1)] = if_else(Response_DF$Response[j] == "Responder", Adjust[91,i],
                                                                       for (k in (2*length(Adjust)):nrow(Response_DF)){
                                                                         Response_DF$Magnitude[i] = if_else(Response_DF$Response[i] == "Neuron", Adjust[110,i], 0.0)}
  )
  )
}
}

} else if (Response_DF$Agonist[i] == "Capsaicin") {
  Response_DF$Magnitude[i] = if_else(Response_DF$Response[i] == "Responder", Adjust[91,i], 0)
} else if (Response_DF$Agonist[i] == "High K") {
  Response_DF$Magnitude[i] = if_else(Response_DF$Response[i] ==  "Neuron", Adjust[110,i], 0)
}
}


if (Response_DF$Agonist[i] == "Capsaicin"){
  