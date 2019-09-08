
install.packages("tidyverse")
library(tidyverse)
install.packages("ggplot2")
library(ggplot2)


##### For non-perfusion agonist additions. Developed by Martha Bhattacharya #####
# import and make data frame and time vectors
well7 = read.csv("M:/Calcium Imaging Analysis/Fluo4_R/well7_fluo4_Results.csv")
# create the time variable; each frame is an image captured 3s apart
t = (well7[,1]*3)-3

# import bkgd files to subtract background
well7bkgd = read.csv("M:/Calcium Imaging Analysis/Fluo4_R/well7bkgd3.csv")
# use three ROIs of background

# average the background fluorescence of each ROI and average the average fluorescence over time
avg_bk = rowMeans(well7bkgd[,2:4])
fluo = well7[,2:ncol(well7)] - avg_bk

# determine the baseline fluorescence of each cell and normalize to their baselines
base = colMeans(fluo[2:6,])
base = data.frame(base)
ratio = data.frame(fluo)
ratio = sweep(fluo, MARGIN=2, FUN = "/", STATS = colMeans(fluo[2:6,]))

plot(t, ratio[,1], xlim = c(0,170), ylim = c(-1,4))



# append the time frame to the normalized fluorescence dataframe
ratio_t = cbind(t,ratio)

# make plots of data
CaImaging = ggplot(data = ratio_t, aes(x = ratio_t$t, y = ratio_t$Mean199)) + geom_point(alpha = 0.5, color = "navy") + coord_cartesian(xlim = c(0,330), ylim = c(0,3)) + labs(x = "Time (Frame*3s)", y = "Normalized Fluorescence (a.u.)", title = "In vitro DRG Calcium Imaging Responses to Various Agonists") + theme(plot.title = element_text(hjust = 0.5)) + theme_light() + geom_hline(yintercept = 1.0, color = "black", linetype = "solid") + geom_vline(xintercept = 27, color = "red", linetype = "dashed") + geom_vline(xintercept = 207, color = "green", linetype = "dashed") + geom_vline(xintercept = 267, color = "purple", linetype = "dashed")

# View the plot
CaImaging

# 10 is a responder. 12? 13? 14? 16?




cellid = colnames(ratio_t)
# extract information from each cell
# first create a list of cell numbers
cellnum = as.numeric(gsub("Mean", "", colnames(ratio)))

peak1 = vector("double", ncol(ratio))
cap1 = vector("double",ncol(ratio))
hik1 = vector("double",ncol(ratio))
for (i in seq_along(ratio)) {
  peak1[i] = max(ratio[11:30,i])
  cap1[i] = max(ratio[71:85,i])
  hik1[i] = max(ratio[91:105,i])
}
allpeaks = data.frame(cbind(cellnum,peak1,cap1,hik1))

##### Data acquisition and prep for perfusion experiment analysis. Developed by Erik Larsen #####


# import and make data frame and time vectors
perf_1 = read.csv("M:/Martha Bhattacharya's Data/calcium imaging/aDRG prep optimization/perf1_table.csv")
perf_1bkg = read.csv("M:/Martha Bhattacharya's Data/calcium imaging/aDRG prep optimization/perf1_bkgd3.csv")


# For each experiment, these variables must be input to perform analysis:

# Frame rate
# Which agonists were used/the order in which they were used
# The time when each agonist was turned on and off (exposure)
# The threshold above normalized baseline which will determine a valid response to agonist exposure


# Define the time variable vector; time was 1.5s / frame, so multiply by 1.5 for time in seconds instead of frames
t = ((perf_1[,1])*1.5)
t = (1:nrow(perf_1))*1.5

# Remove the frame stamp column vector
perf_1 = perf_1[,-1]
perf_1bkg = perf_1bkg[,-1]

# Replace the frame stamp column with the time variable vector (column)
perf_1 = cbind(t, perf_1)
perf_1bkg = cbind(t, perf_1bkg)

# Take the mean of the bkgd for each of three ROIs, then use that mean to normalize the fluorescence in each ROI of the fluorescence data frame
avg_bkgd = rowMeans(perf_1bkg[,2:4])

# (Normalize to background by subtracting the three background ROIs (averaged) from the measured ROIs)
FLUO = perf_1[,2:ncol(perf_1)] - avg_bkgd


# Set the measured cell fluorescences to a value of 1 at time 0 for each cell. This normalizes every frame to the first (1) for all cells/columns
for (i in 1:ncol(FLUO)){
  FLUO[,i] = (FLUO[,i]) / (FLUO[1,i])
}


# Find the minimum, maximum for plotting scales and removing artefacts
min(FLUO[,])
max(FLUO[,])
which(FLUO == min(FLUO), arr.ind = TRUE)
which(FLUO == max(FLUO), arr.ind = TRUE)

# careful with "dev.off()". Avoid when able. COMMENTED OUT
#dev.off()
#plot(t, FLUO[,206], xlim = c(0,331), ylim = c(0,21))




# Create a vector list of agonists in the order they were added in the experiment.

## For 7/12/19, the basis of the development of this script, CQ was added for 15 frames at frame 20, then again from frame 90 to 140; Capsaicin was turned on from frame 140-180; High Potassium from 180-220 ##

# There is a rough, 8s (or ~5 frame) delay for liquid to reach the coverslip from its bottle near the microscope. Once one agonist was turned on, it was the only agonist to which all the cells were exposed.

agonist = c("CQ", "CQ", "Cap", "HiK")

# Create a matrix containing the frame number at which the agonist was added and turned off
exposure_time = cbind(c(25, 95, 145, 185), c(40, 145, 185, 220))

# Rename the start and end columns
colnames(exposure_time) = c("Start", "End")

# Rename the rows with the agonists
rownames(exposure_time) = c(agonist)

# For an example, find the fluorescence recorded over the given interval of one of the agonists in a given cell
FLUO[exposure_time[1,1]:exposure_time[1,2],59]

# Create a data frame to house all the agonists of each experiment along with the percent threshold we deem is worthy of a salient response
Agonists_df = data.frame(c(agonist), c(0.10, 0.10, .05, .05))
colnames(Agonists_df) = c("Agonist", "Percent Response Threshold")


Analysis = function(Agonists_information, Exposure_Information, Fluorescence_Information){


##### Establish Baselines, Thresholds, and Determine Responses for Classification #####

  # These values depend on the onset of agonist in each experiment and subsequent responses
  # Set an arbitrary baseline average FL (pre-agonist) for all cells prior to addition of the first agonist in the 28th frame (added at frame 20 + 8 frame delay)
  baseline_df = matrix(nrow = length(agonist), ncol = length(FLUO))

  # Take the average fluorescence over the first few frames for the first agonist
  baseline_df[1,] = (colSums(FLUO[1:exposure_time[1,1], 1:length(FLUO)], dims = 1))/exposure_time[1,1]

  # Iterate through to repeat for each agonist, except for High Potassium
  for (i in 2:(length(agonist)-1))
    baseline_df[i,] = colSums(FLUO[(exposure_time[i,1]-10):exposure_time[4,1], 1:length(FLUO)], dims = 1)/10

  # Add in the High Potassium threshold
  for (i in length(agonist))
    baseline_df[i,] = colSums(FLUO[(exposure_time[i,1]-2):exposure_time[i,1], 1:length(FLUO)], dims = 1)/2


  # Calculate the response threshold for the first agonist (at frame 25; 5 frames after initiating exposure)
  # Include the thresholds for each agonist in the Agonists_df
  Agonist_threshold_df = matrix(nrow = length(agonist), ncol = length(FLUO))

  # Keep the same cell #
  rownames(Agonist_threshold_df) = c(agonist)

  # Append the threshold values of each cell
  for (i in 1:length(agonist))
    Agonist_threshold_df[i,] = baseline_df[i,]*Agonists_df[i,2] + baseline_df[i,]


  # Determine responders to the first agonist
  # Create an empty vector to fill with whether each cell (ROI) responded to an agonist
  Agonist_responses = data.frame()

  # Fill the vector with cells that responded to the agonist or not, based on the pre-determined, arbitrary threshold of what constitutes a "response"
  for (i in 1:ncol(Agonist_threshold_df)){
    for (j in 1:nrow(Agonist_threshold_df)){
      Agonist_responses[j,i] = if_else(any(FLUO[exposure_time[j,1]:exposure_time[j,2],i] > Agonist_threshold_df[1,i]), "Responder", "Non-responder")
    }
  }


  # Make response calls to High Potassium addition 2 frames after exposure
  for (i in 1:ncol(Agonist_threshold_df)){
    Agonist_responses[4,i] = if_else(Agonist_responses[4,i] == "Responder", "Neuron", 
                                                   if_else(FLUO[182,i] - baseline_df[4,i] > Agonist_threshold_df[4,i], "Neuron", "Non-neuron"))
  }

  # Replace false-positive responses to an agonist
  for (i in 1:length(agonist)) {
    Agonist_responses[i,which(Agonist_responses[i,] == "Responder" & Agonist_responses[length(agonist),] == "Non-neuron")] = "Non-responder"
  }

  # See how many cells are classified as "responders", which ones they are, and the percentages of responses
  for (i in 1:length(agonist)){
    # How many responders?
    print(length(which(Agonist_responses[i,] == "Responder")))
  }

  for (i in 1:length(agonist)){
    # How many neurons?
    print(length(which(Agonist_responses[i,] == "Neuron")))
  }
  
  for (i in 1:length(agonist)){
    # Which are responders?
    print(which(Agonist_responses[i,] == "Responder"))
  }

  for (i in 1:length(agonist)){
    # Which are neurons?
    print(which(Agonist_responses[i,] == "Neuron")) 
  }
  
return(for (i in 1:length(agonist)){
      # Determine percentages
      # Of responders to each agonist; % responders of all neurons, and percent of neurons in the field
      print((length(which(Agonist_responses[i,] == "Responder"))/length(which(Agonist_responses[4,] == "Neuron")))*100)
      print((length(which(Agonist_responses[i,] == "Neuron"))/length(FLUO))*100)
      }
    ) 
}

Analysis(Agonists_information = Agonists_df, Exposure_Information = exposure_time, Fluorescence_Information = FLUO)




which(Agonist_responses[1,] == "Responder" & Agonist_responses[4,] == "Non-neuron")
which(Agonist_responses[2,] == "Responder" & Agonist_responses[4,] == "Non-neuron")
which(Agonist_responses[3,] == "Responder" & Agonist_responses[4,] == "Non-neuron")



# Create a decision data frame
Response_calls = rbind(Agonist_responses[4,])
Response_calls = transpose(Response_calls)



###### *IN PROGRESS* Create a vector variable that gives the quantitative length of time for a High Potassium Response ####
High_K_duration = vector()

for (i in 1:length(Neural_threshold)){
  for (j in 1:length(t)){
    for (k in 182:j) {
      if ((FLUO[(182:j),i]) > FLUO[182,i]) {
        High_K_duration[i] = length(k)}
    }
  }
}
}

High_K_duration[i] = 
  if_else(FLUO[(110:j),i]) > FLUO[110,i], length(t), 
}

}


## IN PROGRESS GARBAGE? ##
dummy = vector()
for (i in 1:length(High_K_responses)){
  if (High_K_responses[i] == "Neuron") {
    dummy[i] = FLUO[35,i]
  } else {
    dummy[i] = 0
  }
}

length(dummy[dummy] != 0)
which(dummy == 0)

which(Agonist_responses[4,] == "Non-neuron")
which(Agonist_responses[1:3,] == "Responder")


##### *IN PROGRESS* Create a new dataframe that will store binary character response classifications of each cell to each reagent along with the response magnitude and duration #####


# Create a dataframe containing the responses of each cell for each agonist





Responses = data.frame()

for (i in 1:length(FLUO)) {
  # 1st Ag Mag
  Responses[i,1] = max(FLUO[20:140,i])
  # 2nd Ag Mag
  Responses[i,2] = max(FLUO[142:179,i])
  # High K Mag
  Responses[i,3] = max(FLUO[180:length(FLUO),i])
  # 1st Ag Duration
  Responses[i,4] = FLUO[(c(11:13)/3),1]i]
# 2nd Ag Duration
Responses[i,5] = FLUO[,i]
# High K Duration
Responses[i,6] = FLUO[,i]
}





Response_DF = data.frame(First_agonist = c(First_agonist_responses), First_ag_Mag = Responses[,1], Capsaicin = c(Capsaicin_responses), Capsaicin_Mag = Responses[,2], High_K = c(High_K_responses), High_K_Mag = Responses[,3])