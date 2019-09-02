#### Larsen Function Library

## Function creating z-scores from a vector of observations, x, the mean of the vector, mx, and the standard deviation of the vector, sd.
zscore = function(x,mx) {return((x-mx)/sd(x))}


## Function that computes the Median of MLB HRs over last 12 Seasons
# Function that finds the position of the median value in a vector
median_fun = function(x) {return(round(length(order(x))/2))}
# Variable for ordering the vector of interest (ML HRs)
ix = order(ML_HR)
# Optional step that calls the order of the vector to make sure it's ordered
ML_HR[ix]
# Creates a variable that stores the values of the vector from lowest to highest
SL = ML_HR[ix]
# Indexes the median value of the vector from the sorted list
SL[median_fun(ML_HR)]

### Note the test function for above (ML_HR)


## Another Function that gives the median of a vector of values
med = function(x) {return(median(x, na.rm = F))}
# Gives the median of ML HR
med(ML_HR)




## Create a normalizing function that takes a vector of numeric values, and for each value, subtracts the minimum value in the vector and divides by the range of values in the vector
normalize = function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


##### Calcium Imaging Functions #####


# To obtain variables to then input into the Data_prep function



# To Analyze Calcium Imaging Data #
Data_prep = function(Adjust, number_agonists, ideal_agonist_responses) {
  # Adjust the baseline fluorescence for each cell to a value of 1 at time 0 in the new response dataframe
  for (i in 1:ncol(ratio)){
    for (j in 1:nrow(ratio)){
      Adjust[j,i] = if_else(ratio[1,i] < 1, ratio[j,i]+(1-ratio[1,i]), ratio[j,i]+(1-ratio[1,i])) 
    }
  }
  
  
  # Account for any artefactual cells with negative values
  for (i in 1:ncol(Adjust)){
    for (j in 1:nrow(Adjust)){
      Adjust[j,i] = if_else(Adjust[j,i] < 0, Adjust[j,i] + 1, Adjust[j,i])
    }
  }
  
  ##### First Agonist #####
  # Set an arbitrary baseline average FL (pre-agonist) for all cells prior to addition of the first agonist in the 9th frame
  # Take the average fluorescence over the first few frames
  ratio_avg_baseline = (colSums(Adjust[1:9, 1:length(Adjust)], dims = 1))/9
  
  
  # Calculate the response threshold for the first agonist (at frame 8; at exposure)
  First_agonist_threshold = ratio_avg_baseline*0.10 + ratio_avg_baseline
  
  
  # Determine responders to the first agonist
  # Create an empty vector to fill with strings, classifying cells as either responders to the given agonist or not
  First_agonist_responses = vector()
  
  # Fill the vector with cells that responded to the agonist or not, based on the pre-determined, arbitrary threshold of what constitutes a "response"
  for (i in 1:length(First_agonist_threshold)){
    First_agonist_responses[i] = if_else(Adjust[25,i] >= First_agonist_threshold[i], "Responder", "Non-responder")
  }
  
  
  Adjust[25,] - ratio_avg_baseline[]
  
  
  # See how many and also which cells are classified as "responders"
  length(First_agonist_responses[First_agonist_responses == "Responder"])
  which(First_agonist_responses == "Responder")
  
  # Percent first agonist responders
  first_agonist_responders = (length(First_agonist_responses[First_agonist_responses == "Responder"]))/ncol(Adjust)*100
  first_agonist_responders
  
  
  First_agonist_duration = vector()
  
  # Establish an average response above which decay in response will be measured as a percentage
  First_agonist_decay_average = (colSums(Adjust[30:44,1:length(Adjust)], dims = 1))/15
  
  # Find an ideal neural response to the first agonist as a model response by cycling through the plots in the new directory on the server
  #Adjust[,14]
  
  # Use that plot number to find the ideal response's maximum response
  max(Adjust[10:30,])
  
  # Find at what frame that response occurs
  which(Adjust[,] == max(Adjust[10:30,]))
  
  # Find the decay (as a percentage from the onset (after peak) of the average)
  ((Adjust[30,] - First_agonist_decay_average[]) / Adjust[30,]*100)
  # In this case, the response decays 14% from peak response to the average or baseline following addition of High Potassium
  
  
  # Find the peak of decay (as a percentage from the onset of the average)
  ((Adjust[30,] / Adjust[35,])-1)*100
  Adjust[35:44,] > First_agonist_decay_average[]
  # In this case, the response decays 14.9% from peak response to the average or baseline following addition of first agonist

  
  
  
  # Find the duration of responses to the first agonist; how long does the response last?
  
  for (i in 1:length(First_agonist_responses)){
    if (First_agonist_responses[i] == "Responder") {
      First_agonist_duration[i] = length(which(Adjust[10:35,i] > First_agonist_decay_average[i]))
    } else {
      First_agonist_duration[i] = 0
    }
  }
  
  
  # Find how many Responders to the first agonist had decays
  length(First_agonist_duration[First_agonist_duration > 0])
  # Find which Responders to the first agonist had decays
  which(First_agonist_duration > 0) & which(First_agonist_responses == "Responder")
  
  
  
  First_agonist_responses[]
  First_agonist_duration[]
  
  Adjust[25,] / Adjust[30,]
  
  
  
  ##### Second Agonist (CQ or Histamine) #####
  
  
  # Determine responders to capsaicin
  # Establish a baseline, accounting for the mean fluorescence in the 9 frames leading up to the second agonist addition
  Second_agonist_baseline = (colSums(Adjust[36:44,1:length(Adjust)], dims = 1))/9
  
  
  # Set the threshold 1 frame after exposure
  Second_agonist_threshold = Second_agonist_baseline*0.10 + Second_agonist_baseline
  
  
  # Create an empty vector to fill with strings, classifying cells as either responders to the given agonist or not
  Second_agonist_responses = vector()
  
  
  
  # Make response calls after addition in the 45th frame
  for (i in 1:length(Second_agonist_threshold)){
    Second_agonist_responses[i] = if_else(Adjust[45,i] > Second_agonist_threshold[i], "Responder", "Non-responder")
  }
  
  # See how many and also which cells are classified as "responders"
  length(Second_agonist_responses[Second_agonist_responses == "Responder"])
  which(Second_agonist_responses == "Responder")
  
  second_agonist_responders = (length(Second_agonist_responses[Second_agonist_responses == "Responder"]))/ncol(Adjust)*100
  second_agonist_responders
  
  
  
  Second_agonist_duration = vector()
  
  # Establish an average response where the rise in response asympotitcally slows (approaches ceiling); will be measured as a percentage
  Second_agonist_decay_average = (colSums(Adjust[65:79,1:length(Adjust)], dims = 1))/15
  
  # Find an ideal response to the Second agonist as a model response
  #Adjust[,147]
  
  
  # Find the ideal response's maximum response
  max(Adjust[65:79,])
  
  # Find at what frame that response occurs
  which(Adjust[,] == max(Adjust[65:79,]))
  
  # Find the decay (as a percentage from the onset (after peak) of the average)
  ((Adjust[66,] - Second_agonist_decay_average[]) / Adjust[66,]*100)
  # In this case, the response decays 14% from peak response to the average or baseline following addition of High Potassium
  
  
  # Find the peak of decay (as a percentage from the onset of the average)
  ((Adjust[66,] / Adjust[66,])-1)*100
  Adjust[60:70,] > Second_agonist_decay_average[]
  # In this case, the response decays 14.9% from peak response to the average or baseline following addition of second agonist
  
  
  
  # Find the duration of responses to the second agonist
  
  for (i in 1:length(Second_agonist_responses)){
    if (Second_agonist_responses[i] == "Responder") {
      Second_agonist_duration[i] = length(which(Adjust[60:70,i] > Second_agonist_decay_average[i]))
    } else {
      Second_agonist_duration[i] = 0 
    }
  }
  
  
  
  ##### Capsaiscin Threshold for Classification #####
  
  
  # Determine responders to capsaicin
  # Establish a baseline, accounting for the mean fluorescence in the 8 frames leading up to the second agonist addition
  Capsaicin_baseline = (colSums(Adjust[71:79,1:length(Adjust)], dims = 1))/9
  
  
  # Set the threshold at frame 71; 1 frame after exposure
  Capsaicin_agonist_threshold = Capsaicin_agonist_baseline*0.10 + Capsaicin_agonist_baseline
  
  
  # Create an empty vector to fill with strings, classifying cells as either responders to the given agonist or not
  Capsaicin_agonist_responses = vector()
  
  
  # Make response calls after addition in the 80th frame
  for (i in 1:length(Capsaicin_agonist_threshold)){
    Capsaicin_agonist_responses[i] = if_else(Adjust[90,i] > Capsaicin_agonist_threshold[i], "Responder", "Non-responder")
  }
  
  # See how many and also which cells are classified as "responders"
  length(Capsaicin_agonist_responses[Second_agonist_responses == "Responder"])
  which(Capsaicin_agonist_responses == "Responder")
  
  Capsaicin_agonist_responders = (length(Capsaicin_agonist_responses[Capsaicin_agonist_responses == "Responder"]))/ncol(Adjust)*100
  Capsaicin_agonist_responders
  
  
  
  Capsaicin_agonist_duration = vector()
  
  # Establish an average response where the rise in response asympotitcally slows (approaches ceiling); will be measured as a percentage
  Capsaicin_agonist_decay_average = (colSums(Adjust[90:95,1:length(Adjust)], dims = 1))/5
  
  # Find an ideal response to the Second agonist as a model response
  #Adjust[,147]
  
  
  # Find the ideal response's maximum response
  max(Adjust[85:95,])
  
  # Find at what frame that response occurs
  which(Adjust[,] == max(Adjust[85:95,]))
  
  # Find the decay (as a percentage from the onset (after peak) of the average)
  ((Adjust[90,] - Second_agonist_decay_average[]) / Adjust[90,]*100)
  # In this case, the response decays 14% from peak response to the average or baseline following addition of High Potassium
  
  
  # Find the peak of decay (as a percentage from the onset of the average)
  ((Adjust[90,] / Adjust[90,])-1)*100
  Adjust[90:95,] > Second_agonist_decay_average[]
  # In this case, the response decays 14.9% from peak response to the average or baseline following addition of second agonist
  
  
  
  # Find the duration of responses to Capsaicin
  
  for (i in 1:length(Capsaicin_agonist_responses)){
    if (Capsaicin_agonist_responses[i] == "Responder") {
      Capsaicin_agonist_duration[i] = length(which(Adjust[85:104,i] > Capsaicin_agonist_decay_average[i]))
    } else {
      Capsaicin_agonist_duration[i] = 0 
    }
  }
  
  
  
  
  ##### Neural Threshold for Classification #####
  
  
  # Determine neurons
  # Establish a baseline, accounting for the mean in the 8 frames leading up to the High K addition
  Neural_baseline = (colSums(Adjust[80:86,1:length(Adjust)], dims = 1))/7
  
  
  # Set the threshold at 5% above the baseline established 8 frames before High K addition
  Neural_threshold = Neural_baseline*0.05 + Adjust[88,1:length(Adjust)]
  
  
  # Create an empty vector to fill with strings, classifying cells as either neurons or not
  High_K_responses = vector()
  
  # Make response calls in the 110th frame to High K addition in the 109th frame
  for (i in 1:length(Neural_threshold)){
    High_K_responses[i] = if_else(Second_agonist_responses[i] == "Responder", "Neuron", 
                                  if_else(Adjust[89,i] - Neural_baseline[i] > 0.5, "Neuron", "Non-neuron"))
  }
  
  
  # See how many and also which cells are classified as "responders"
  length(High_K_responses[High_K_responses == "Responder"])
  which(High_K_responses == "Responder")
  
  # What percent of the ROIs are neurons?
  high_K_responders = (length(High_K_responses[High_K_responses == "Neuron"]))/ncol(Adjust)*100
  high_K_responders
  
  
  
  
  High_K_duration = vector()
  
  # Establish an average response above which decay will be measured as a percentage in the first half of the time after addition (to total time)
  High_K_decay_average = (colSums(Adjust[92:103,1:length(Adjust)], dims = 1))/12
  
  # Find an ideal neural response to High Potassium as a model response
  #Adjust[,7]
  
  # Find the ideal response's maximum response
  max(Adjust[,7])
  
  # Find at what frame that response occurs
  which(Adjust[,7] == max(Adjust[,7]))
  
  # Find the decay (as a percentage from the onset (after peak) of the average)
  ((Adjust[92,7] - High_K_decay_average[7]) / Adjust[92,7]*100)
  # In this case, the response decays 12% from peak response to the average or baseline following addition of High Potassium
  
  
  
  # Find the duration of neural responses to High Potassium
  
  for (i in 1:length(High_K_responses)){
    if (High_K_responses[i] == "Neuron") {
      High_K_duration[i] = length(which(Adjust[90:103,i] > High_K_decay_average[i]))
    } else {
      High_K_duration[i] = 0 }
  }
  
  
  
  
  
  # See how many and also which cells are classified as "Neurons"
  length(High_K_responses[High_K_responses == "Neuron"])
  which(High_K_responses == "Neuron")
  
  
  # How many cells "respond" to the first agonist, yet aren't neurons?
  length(which(High_K_responses == "Non-neuron" & First_agonist_responses == "Responder"))
  
  
  ##### Create a new dataframe that will store binary character response classifications of each cell to each reagent along with the response magnitude and duration #####
  
  
  # Create a dataframe containing the responses of each cell for each agonist
  
  
  
  
  
  Responses = data.frame()
  
  for (i in 1:length(Adjust)) {
    # 1st Ag Mag
    Responses[i,1] = Adjust[9,i]
    # 2nd Ag Mag
    Responses[i,2] = Adjust [69,i]
    # High K Mag
    Responses[i,3] = Adjust [89,i]
    # 1st Ag Duration
    Responses[i,4] = First_agonist_duration[i]
    # 2nd Ag Duration
    Responses[i,5] = Second_agonist_duration[i]
    # High K Duration
    Responses[i,6] = High_K_duration[i]
    # 1st Ag Class.
    Responses[i,7] = First_agonist_responses[i]
    # 2nd Ag Class.
    Responses[i,8] = Second_agonist_responses[i]
    # High K Class.
    Responses[i,9] = High_K_responses[i]
  }
  
  colnames(Responses) = c("1st Agonist Response Magnitude", "2nd Agonist Response Magnitude", "High K Response Magnitude",
                          "1st Agonist Response Duration", "2nd Agonist Response Duration", "High K Response Duration",
                          "1st Agonist Responder", "2nd Agonist Responder", "High K Responder")
  
  
  
  
  
  Response_DF = data.frame(First_agonist = c(First_agonist_responses), First_ag_Mag = Responses[,1], Capsaicin = c(Capsaicin_responses), Capsaicin_Mag = Responses[,2], High_K = c(High_K_responses), High_K_Mag = Responses[,3])
  
  
  
  
  
  
  
  
  return (Responses, Agonist_proportions)
}
