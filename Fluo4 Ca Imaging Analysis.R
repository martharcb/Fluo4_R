

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
# Number/which agonists were used
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

# careful with "dev.off()". Avoid when able.
#dev.off()
#plot(t, FLUO[,206], xlim = c(0,331), ylim = c(0,21))




# Create a vector list of agonists in the order they were added in the experiment.

# For 7/12/19, the basis of the development of this script, CQ was added for 15 frames at frame 20, then again from frame 90 to 140; Capsaicin was turned on from frame 140-180; High Potassium from 180-220

# There is a rough, 8s or 8 frame delay for liquid to reach the coverslip from its bottle near the microscope. Once one agonist was turned on, it was the only agonist to which all the cells were exposed.

agonist = c("CQ", "CQ", "Cap", "HiK")

# Create a matrix containing the frame number at which the agonist was added and turned off
exposure_time = cbind(c(28, 98, 148, 188), c(43, 148, 188, 196))

# Rename the start and end columns
colnames(exposure_time) = c("Start", "End")

# Rename the rows with the agonists
rownames(exposure_time) = c(agonist)

# For an example, find the fluorescence recorded over the given interval of one of the agonists in a given cell
FLUO[exposure_time[1,1]:exposure_time[1,2],59]

# Create a data frame to house all the agonists of each experiment along with the percent threshold we deem is worthy of a salient response
Agonists_df = data.frame(c(agonist), c(0.10, 0.10, .05, .05))
colnames(Agonists_df) = c("Agonist", "Percent Response Threshold")





##### First Agonist Threshold for Classification #####

# These values depend on the onset of agonist in each experiment and subsequent responses


# Set an arbitrary baseline average FL (pre-agonist) for all cells prior to addition of the first agonist in the 28th frame (added at frame 20 + 8 frame delay)
# Take the average fluorescence over the first few frames
baseline_df = matrix(nrow = length(agonist), ncol = length(FLUO))
  
baseline_df[1,] = (colSums(FLUO[1:exposure_time[1,1], 1:length(FLUO)], dims = 1))/exposure_time[1,1]


# Calculate the response threshold for the first agonist (at frame 28; 8 frames after exposure)

# Include the thresholds for each agonist in the Agonists_df
Agonist_threshold_df = matrix(nrow = length(agonist), ncol = length(FLUO))

# Keep the same cell #
rownames(Agonist_threshold_df) = c(agonist)

# Append the threshold values of each cell
Agonist_threshold_df[1,] = baseline_df[1,]*Agonists_df[1,2] + baseline_df[1,]


# Determine responders to the first agonist
# Create an empty vector to fill with whether each cell (ROI) responded to an agonist
Agonist_responses = data.frame()


# Fill the vector with cells that responded to the agonist or not, based on the pre-determined, arbitrary threshold of what constitutes a "response"
for (i in 1:ncol(Agonist_threshold_df)){
  Agonist_responses[1,i] = if_else(any(FLUO[exposure_time[1,1]:exposure_time[1,2],i] > Agonist_threshold_df[1,i]), "Responder", "Non-responder")
}


# For checking
#FLUO[11,38] - baseline[38]


# See how many (and also which cells) are classified as "responders"
length(Agonist_responses[1][Agonist_responses == "Responder"])
which(Agonist_responses == "Responder")




##### Second Agonist Threshold for Classification #####

# Determine responders to the second addition of CQ
baseline_df[2,] = (colSums(FLUO[(exposure_time[2,1]-10):exposure_time[2,1], 1:length(FLUO)], dims = 1))/10


# Calculate the response threshold for the second agonist (at frame 98; 8 frames after exposure)

# Append the threshold values of each cell
Agonist_threshold_df[2,] = baseline_df[2,]*Agonists_df[2,2] + baseline_df[2,]


# Fill the vector with cells that responded to the agonist or not, based on the pre-determined, arbitrary threshold of what constitutes a "response"
for (i in 1:ncol(Agonist_threshold_df)){
  Agonist_responses[2,i] = if_else(any(FLUO[exposure_time[2,1]:exposure_time[2,2],i] > Agonist_threshold_df[2,i]), "Responder", "Non-responder")
}


# For checking
#FLUO[11,38] - baseline[38]


# See how many (and also which cells) are classified as "responders"
length(Agonist_responses[2][Agonist_responses == "Responder"])
which(Agonist_responses == "Responder")




##### Capsaiscin Threshold for Classification #####



# Determine responders to capsaicin
# Establish a baseline, accounting for the mean fluorescence in the 10 frames leading up to the Capsaicin addition
baseline_df[3,] = (colSums(FLUO[130:140,1:length(FLUO)], dims = 1))/10


# Set the threshold at frame 142; 2 frames after exposure
Agonist_threshold_df[3,] = baseline_df[3,]*Agonists_df[3,2] + FLUO[142,1:length(FLUO)]


# Make response calls after addition in the 142nd frame
for (i in 1:length(Agonist_threshold_df)){
  Agonist_responses[3,i] = if_else(any(FLUO[142:179,i] - baseline_df[3,i] > Agonists_df[3,2]), "Responder", "Non-responder")
}

# See how many and also which cells are classified as "responders"
length(Agonist_responses[3][Agonist_responses == "Responder"])
which(Agonist_responses == "Responder")




##### Neural Threshold for Classification #####


# Determine neurons
# Establish a baseline, accounting for the mean in the 8 frames leading up to the High K addition
baseline_df[4,] = (colSums(FLUO[172:180,1:length(FLUO)], dims = 1))/8


# Set the threshold at 5% above the baseline established 2 frames before High K addition
Agonist_threshold_df[length(agonist),i] = baseline_df[length(agonist),]*Agonists_df[4,2] + FLUO[178,1:length(FLUO)]


# Make response calls to High Potassium addition 2 frames after exposure
for (i in 1:length(Agonist_threshold_df)){
  Agonist_responses[length(agonist),i] = if_else(Agonist_responses[3,i] == "Responder", "Neuron", 
                                if_else(FLUO[182,i] - baseline_df[length(agonist),i] > Agonists_df[4,2], "Neuron", "Non-neuron"))
}

# *IN PROGRESS* Create a vector variable that gives the quantitative length of time for a High Potassium Response
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

# See how many and also which cells are classified as "Neurons"    *IN PROGRESS*
length(which(Agonist_responses == "Neuron"))
which(Agonist_responses == "Neuron")


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

which(Agonist_responses == "Non-neuron" & Agonist_responses == "Responder")


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

















# Make individual plots as separate png files; these are response profiles for one cell (ROI), each. These give rough ideas of responses by each ROI
dir.create("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fluo4 Images\\perf1plots")
for (i in seq_along(FLUO[,ncol(FLUO)])) {
  #tempname = paste("cell_",i,".png", sep="")
  outpath = file.path("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fluo4 Images\\perf1plots", paste("cell_",i,".png", sep=""))
  dev.off()
  png(filename = outpath)
  plot(t, FLUO[,i], xlim = c(0,331), ylim = c(0,6.1))
}

# Add the time vector again
FLUO = cbind(t, FLUO)

# Make plots of perfusion data
CaImaging = ggplot(data = FLUO, aes(x = FLUO$t, y = FLUO$Mean70)) + geom_vline(xintercept = 30, color = "forest green") + geom_vline(xintercept = 52.5, color = "forest green") + geom_vline(xintercept = 135, color = "forest green") + geom_vline(xintercept = 208.5, color = "forest green") + geom_vline(xintercept = 210, color = "red") + geom_vline(xintercept = 238.5, color = "red") + geom_vline(xintercept = 285, color = "purple") + geom_point(alpha = 0.4, color = "blue") + coord_cartesian(xlim = c(0,331), ylim = c(0,6.1)) + labs(x = "Time (Frame; 1.5s/Frame)", y = "Normalized Fluorescence (a.u.)", title = "In vitro WT DRG Calcium Imaging Responses to CQ, CQ, Cap") + theme(plot.title = element_text(hjust = 0.5)) + theme_light() + geom_hline(yintercept = 1.0, color = "black", linetype = "solid") 

dev.off()
CaImaging




#+ geom_rect(aes(xmin = 30, ymin = -Inf, xmax = 52.5, ymax = Inf), fill = "forest green", color = "forest green", size = 1) + geom_rect(aes(xmin = 135, ymin = -Inf, xmax = 210, ymax = Inf), fill = "forest green", color = "forest green", size = 1) + geom_rect(aes(xmin = 210, ymin = -Inf, xmax = 285, ymax = Inf), fill = "red", color = "red", size = 1) + geom_rect(aes(xmin = 285, ymin = -Inf, xmax = 330, ymax = Inf), fill = "purple", color = "purple", size = 1)































##### Data acquisition and prep for perfusion experiment 2 analysis #####


# import and make data frame and time vectors
perf_2 = read.csv("M:/Martha Bhattacharya's Data/calcium imaging/aDRG prep optimization/2019-07-11/tiff export/perf2 TIFF export/perf2_table.csv")
perf_2bkg = read.csv("M:/Martha Bhattacharya's Data/calcium imaging/aDRG prep optimization/perf2_bkgd3.csv")


# Define the time variable vector
t = ((perf_2[,1])*1.5)
t = (1:nrow(perf_2))*1.5

# Remove the frame stamp vector
perf_2 = perf_2[,-1]
perf_2bkg = perf_2bkg[,-1]

# Replace the frame stamp vector with the time variable vector
perf_2 = cbind(t, perf_2)
perf_2bkg = cbind(t, perf_2bkg)

# Take the mean of the bkgd for each of three ROIs, then use to normalize the fluorescence in each ROI of signal data frame
avg_bkgd = rowMeans(perf_2bkg[,2:4])



FLUO = perf_2[,2:ncol(perf_2)] - avg_bkgd

# Normalize the fluorescences to 1 at time 0 for all frames for each cell
for (i in 1:ncol(FLUO)){
  FLUO[,i] = (FLUO[,i]) / (FLUO[1,i])
}


# Find the minimum, maximum for plotting scales and removing artefacts
min(FLUO[,])
max(FLUO[,])
which(FLUO == min(FLUO), arr.ind = TRUE)
which(FLUO == max(FLUO), arr.ind = TRUE)

dev.off()
plot(t, FLUO[,206], xlim = c(0,331), ylim = c(0,21))




# Make individual plots as separate png files; these are response profiles for one cell (ROI), each. These give rough ideas of responses by each ROI
dir.create("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fluo4 Images\\perf2plots")
for (i in seq_along(FLUO[,ncol(FLUO)])) {
  #tempname = paste("cell_",i,".png", sep="")
  outpath = file.path("M:\\Erik\\Data\\Calcium Imaging Analysis\\Fluo4 Images\\perf2plots", paste("cell_",i,".png", sep=""))
  dev.off()
  png(filename = outpath)
  plot(t, FLUO[,i], xlim = c(0,331), ylim = c(0,6))
}

# Add the time vector again
FLUO = cbind(t, FLUO)

# Make plots of perfusion data
CaImaging = ggplot(data = FLUO, aes(x = FLUO$t, y = FLUO$Mean70)) + geom_vline(xintercept = 30, color = "forest green") + geom_vline(xintercept = 52.5, color = "forest green") + geom_vline(xintercept = 135, color = "forest green") + geom_vline(xintercept = 208.5, color = "forest green") + geom_vline(xintercept = 210, color = "red") + geom_vline(xintercept = 238.5, color = "red") + geom_vline(xintercept = 285, color = "purple") + geom_point(alpha = 0.4, color = "blue") + coord_cartesian(xlim = c(0,331), ylim = c(0,6)) + labs(x = "Time (Frame; 1.5s/Frame)", y = "Normalized Fluorescence (a.u.)", title = "In vitro WT DRG Calcium Imaging Responses to CQ, CQ, Cap") + theme(plot.title = element_text(hjust = 0.5)) + theme_light() + geom_hline(yintercept = 1.0, color = "black", linetype = "solid") 

dev.off()
CaImaging