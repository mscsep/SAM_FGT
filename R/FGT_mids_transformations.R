# Functions for the transformation and manipulation of mids objects in the analysis of the data from the Fear Generalization Task (FGT) in the SAM study
# written by Milou Sep

# Info on how to perform manipulations/calculations on midsobject: 
# https://stackoverflow.com/questions/26667162/perform-operation-on-each-imputed-dataset-in-rs-mice

library(dplyr)

# Function to calculate Means on subset mids object data (bv early trials only..)
means_EML_2way <- function(data_subset){
  data_subset.means<-aggregate(data_subset$FPS,
                               by=list(data_subset$.imp, data_subset$pp, data_subset$Condition, data_subset$trialtype), 
                               FUN=mean)
  names(data_subset.means)<-c(".imp", "pp", "Condition", "trialtype", "FPS")
  # Some checks
  print(unique(data_subset.means$Condition))
  print(unique(data_subset.means$pp))
  print(unique(data_subset.means$.imp))
  
  data_subset.means.mids <- as.mids(data_subset.means) # Make mids.
  return(data_subset.means.mids)
}

# Function to log-transform FPS data if necessary (checked in Assumption Checks)
log.transform.mids <- function(data_set){
  data_set.long<-complete(data_set, action = "long", include = T) # Change mids to long format
  data_set.long$FPS <- log(data_set.long$FPS + 1) # Add 1 to deal with 0-responses
  log.mids.dataset <- as.mids(data_set.long)
  return(log.mids.dataset)
}


# Select one imputed dataset for testing (Note function returns dataframe!)
only_one_imputation <- function(data_mids, imp){
  # 1) Mids to dataframe
  data_mids.long <- complete(data_mids, action = "long", include = T)
  tibble(data_mids.long)
  # 2) Perform manipulations: 
  data_mids.long %>% filter(.imp == imp) ->x
  return(x)
}


# For sensitivity analyses ------------------------------------------------

# Sensitivity analyses to check the effect of potential influential participants:
remove.influential.points <- function(data_mids, influentialpoints){
  # Note numbers should be assigned to influential points e.g. c(1,2,3)
  # Mids to dataframe
  data_mids.long <- complete(data_mids, action = "long", include = T)
  # Perform manipulations: search & exclude rows from influential points (or participants)
  indices <- which(data_mids.long$pp %in% influentialpoints)
  data_mids.long.noinfl <- data_mids.long[-c(indices),]
  # Back to mids.
  data_mids.noinfl <- as.mids(data_mids.long.noinfl)
  return(data_mids.noinfl)
}

# Function to change factor Trial and/or Trialtype to continuous variables 
Within.continuous.mids <- function(data_set){ 
  # This was recommended for analyses of fear gradients by https://www.sciencedirect.com/science/article/pii/S0005791618300612 & http://www.frontiersin.org/Quantitative_Psychology_and_Measurement/10.3389/fpsyg.2015.00652/abstract
  data_set.long<-complete(data_set, action = "long", include = T) # Change mids to long format
  # manipulations:
  data_set.long$trialtype <- as.numeric(data_set.long$trialtype)
  data_set.long$trial <- as.numeric(data_set.long$trial)
  str(data_set.long)
  # transform back to mids
  mids.dataset <- as.mids(data_set.long)
  return(mids.dataset)
}

# Sensitivity analyses on participants that detected the threat stimulus correctly (according to self-report)
Sensitivity_contingency <- function(midsobject){
  long_df<-complete(midsobject, action = "long", include = T)
  pp_to_include <- which(long_df$pp %in% pp_Sensitivity$X1)
  long_df_included <- long_df[c(pp_to_include),]
  mids_included <- as.mids(long_df_included)
  return(mids_included)
}