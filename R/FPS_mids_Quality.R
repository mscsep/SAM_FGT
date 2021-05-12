# Create MIDS-objects from the imputed startle data (EMG) of the Fear Generalization Task (FGT) in the SAM study, and check the imputation quality.
# Written by Milou Sep

# This script can be used to create MIDS objects (a datatype form the mice package) from the output of the FGT imputation script: "FPS_imputation.R"
# The central function in this script is "MAKE_MIDS_Imputed_FGT_EMG_SAM", in which the data is sorted, factors are created and mids objects is formed
# Note This script works with both output of sorttype 'mean' and 'trials' and always creates separate mids objects for day1 and day2.

# For information on the imputation quality checks: https://www.kaggle.com/captcalculator/imputing-missing-data-with-the-mice-package-in-r.

# * Load required packages --------------------------------------------------
# install.packages("mice"); install.packages("car"); utils::install.packages("miceadds", dependencies = TRUE); install.packages("reshape2")
library(mice)
library(car)
library(miceadds)
library(reshape2)
library(stringr)  # Split variable name in multiple columns (info from: https://stackoverflow.com/questions/4350440/split-data-frame-string-column-into-multiple-columns)

# * Load required files -----------------------------------------------------
# Load raw data (umag):
load(file = 'processed_data/INPUT_IMPUTATIE_FGT_Umag.rda', envir = .GlobalEnv, verbose = FALSE) 
# Load imputed data umag (Note M=100 and MAXIT=50):
load(file = 'processed_data/OUPUT_IMPUTATIE_FGT_out_M50_MAXIT100_Umag_Trials.rda', envir = .GlobalEnv, verbose = FALSE) 
load(file= "processed_data/OUPUT_IMPUTATIE_FGT_out_M50_MAXIT100_Umag_AllMeans.rda", envir = .GlobalEnv, verbose = FALSE)

# Check amount of missing data --------------------------------------------
summary(Umag)
n.missing<-sum(is.na(Umag))
# total number of observations
n.observations<-nrow(Umag)*(ncol(Umag)-2) # Note: -2 to remove subject en condition variables from the count
# percent.missing
(n.missing/n.observations)*100  # = 10.39886

# Number of missing values in the imputed trials datasets
(sum(is.na(out_Umag_trials$Trials_imputed[[1]]) ) / n.observations) *100 # = 7,13141   
summary(out_Umag_trials$Trials_imputed[[1]])

# Number of missing values in the imputed means (based on all trials) datasets
summary(out_Umag_all$Means_Completecases) # Note data with missing values
aantal_missings_m.all<-sum(is.na(out_Umag_all$Means_Completecases))
observations_m.all <- nrow(out_Umag_all$Means_Completecases) * ncol(out_Umag_all$Means_Completecases)
(aantal_missings_m.all/observations_m.all) * 100 # = 24.13273%
# To check
summary(out_Umag_all$MergedMeans_imputed[[1]]) # imputed data, no missing values


# Definition of Function Make Mids" --------------------------------------
MAKE_MIDS_Imputed_FPS <- function(input_data, output_data, M, var_pattern, within.factors, number.trialtypes, iti.mean){
  # This function can be used to select FPS responses to cue or context stimuli in the acquisition and test phase of the FGT task
  # Selected data is changed to long format and factors trialtype and trial will be added to the dataset based on original variable names
  # Note 1: number.trialtypes is 2 if only threat and safe trials are present and 3 if threat, safe and new trials are present.
  # Note 2: M is the number of imputations
  
  # Make data.frame with (Subject Name) & (experimental) Condition 
  Subject_Condition <- data.frame(factor(input_data$subject), 
                                  factor(input_data$Condition, levels=c("1","2","3"), 
                                         labels=c("Delayed Stress", "Direct Stress", "No Stress")))
  colnames(Subject_Condition) <- c("pp", "Condition")
  
  # To check if condition row's are identical (if they are both 'factor')
  # identical(factor(Subject_Condition$Smag.Condition), factor(Subject_Condition$output_data..1....m...Condition_i)) # Yes, so both can be used in later dataset
  
  Data_List<-list()  # Define empty list (for main mids object)
  
  for(m in 0:M){ # Loop over m is 1 t/m 101 (including the complete cases dataset, the original data, and 100 imputed datasets)
    
    # Select either Complete Case data (if m is 0) ...
    if (m==0){ Data <- output_data[[2]][ ,grep(pattern=var_pattern, x=names(data.frame(output_data[[2]])))] # Note data.frame() required for names().
    # ... or a imputed datasets (if m is NOT 0).
    } else if (m!=0) {Data <- output_data[[1]][[m]][ ,grep(pattern=var_pattern, x=names(output_data[[1]][[m]]))]}
    
    # Bind data with Subject name & Experimental condition...
    Data <- cbind(Subject_Condition, Data) 
    # ... and change data to long format
    Data_Long <- melt(Data, id.vars = c("pp", "Condition"), value.name = 'FPS')
    
    
  # If the individual trials are analyzed: iti.mean = 0 (so no factors need to be created)
    if(iti.mean == 0){
      
      # Split the variable name in 5 elements on "_", add the 4th element that is created (e.g. "T1" or "S4") as variable "TrialCode" to `datalong`...
      variable.codes<-data.frame(stringr::str_split_fixed(string=Data_Long$variable, pattern="_",n=5)) 
      ## To check results of 'str_split_fixed' for cue acquisition trials:
      # unique(variable.codes[,1]) # sA     [ for cue acq all & EL: Mean] 
      # unique(variable.codes[,2]) # cue    [ for cue acq all & EL: sA] 
      # unique(variable.codes[,3]) # Umag   [ for cue acq all & EL: cue]
      # unique(variable.codes[,4]) # T1, T2, T3, ..etc  S1, S2, S3 etc.    [ for cue acq all & EL: T & S, for ITI E&L is this empty]
      # unique(variable.codes[,5]) # i      [for cue acq all: merged, for cue acq EL: contains E&L, for ITI EL this is E&L}
      
      #  ... and trialtype too:
      if("trialtype" %in% within.factors){
        # If both trialtype and trial are BOTH within factors, which is only the case in datasets for analyses
        # of the acquisition and test phase of cued and contextual FPS, variable.code 4 needs to be divided in letters & digits
        # (for info see https://stackoverflow.com/questions/9756360/split-character-data-into-numbers-and-letters)
        # Numbers (recognized by pattern: "[[:digit:]]") are replaced by nothing (""), as a result only letters remain.
        Data_Long$trialtype<-factor(gsub("[[:digit:]]","", variable.codes[,4])) # Note: gsub(pattern, replacement, x)
        # Idem, only here are 'not numbers' (so the letters) replaced by nothing ("").
        Data_Long$trial<-factor(gsub("[^[:digit:]]","", variable.codes[,4]))
        
        # If trialtype is a within factor, recode trialtypes Threat, Safe & New
        # Recode trialtype & make factor
        if (number.trialtypes == 2){
          # Recode T en S to 0 en 1
          Data_Long$trialtype<-recode(var=Data_Long$trialtype, recodes="'T'= 0; 'S'= 1") 
          # Make trialtype factor.
          Data_Long$trialtype <- factor(Data_Long$trialtype, levels = c("0","1"), labels = c("Threat", "Safe"))
        } else if (number.trialtypes == 3) {
          # Recode T en S en N to 0 en 1 en 2
          Data_Long$trialtype<-recode(var=Data_Long$trialtype, recodes="'T'= 0; 'S'= 1; 'N'=2") 
          # Make trialtype factor.
          Data_Long$trialtype <- factor(Data_Long$trialtype, levels = c("0","1", "2"), labels = c("Threat", "Safe", "New"))
        }
        
        # ... and if trialtype is NO within factor:
      } else if (!("trialtype" %in% within.factors)){
        # If trialtype is not a within factor (and no means are analyzed), which is the case in analyses of HAB & ITI trials:
        Data_Long$trial<-factor(variable.codes[,4])
      }
      
    }
    
    # To check the results of the data manipulations
    # str(Data_Long)
    # unique(Data_Long$trialtype)
    
    Data_List[[m+1]] <- Data_Long # Note m+1 is needed because M starts at 0 in this loop, but 0 does not indicate a position in the list.
  }
  
  # Note do.call applies rbind() to on all list elements (M+1)
  Data_Frame <- do.call(rbind,Data_List)
  # Add column's .imp en .id (both necessary to make mids object, via as.mids()!
  # Add code .imp (required to indicate the number of imputed data set, 0 = with missing values, 1 = imputed dataset 1 etc. (M))
  Data_Frame$.imp <- factor(c(rep(0:M,each=c(nrow(Data_List[[1]])))))
  # Add .id column (that indicates the row in the dataset)
  Data_Frame$.id <- factor(c(1:nrow(Data_Frame)))
  
  # Create mids object(s) for overall dataframe..
  mids.data<-as.mids(long=Data_Frame, .imp=".imp", .id=".id")
  
  # Define Function output    
  # return(Data_List)  # list with M+1 elements (dataframes), the first is the original data (with missing values)
  # return(Data_Frame) # rows from the M+1 dataframes below each other in one dataframe
  return(mids.data) # mids object of variable "Data_Frame"
}

# Call script above, per trial type:

# Acquisition Cued FPS ----------------------------------------------------
Acquisition_Cued_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS( 
  input_data=Umag,output_data = out_Umag_trials, M=100, var_pattern = 'sA_Cue', 
  within.factors=c("trial", "trialtype"), number.trialtypes=2, iti.mean=0)
# Check quality imputations
print(densityplot(Acquisition_Cued_FPS_mids_Trials, ~FPS))
bwplot(Acquisition_Cued_FPS_mids_Trials)
# MS & RG: OK

# Acquisition Context FPS -------------------------------------------------
Acquisition_Context_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag, output_data = out_Umag_trials, M=100, var_pattern = 'sA_Ctx', 
  within.factors=c("trial", "trialtype"), number.trialtypes=2, iti.mean=0)
print(densityplot(Acquisition_Context_FPS_mids_Trials, ~FPS))
bwplot(Acquisition_Context_FPS_mids_Trials)
# MS & RG: OK

# Context Habituation -----------------------------------------------------
HABctx_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag, output_data = out_Umag_trials, M=100, var_pattern = 'sA_HABctx', 
  within.factors=c("trial", "trialtype"), number.trialtypes=3, iti.mean=0)
# No imputations for this trialtype.
# print(densityplot(HABctx_FPS_mids_Trials, ~FPS)) 
# bwplot(HABctx_FPS_mids_Trials)
summary(complete(HABctx_FPS_mids_Trials)) # 72 missings (12 participants, 3 categories, 3trials per category)
# 12*2*3
summary(Umag) # for 12 participants are all habctx trials missing (= more then 2/3 missing), so no imputations were performed.

# Test Cued FPS -----------------------------------------------------------
Test_Cued_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag, output_data = out_Umag_trials, M=100, var_pattern = 'sG_Cue', 
  within.factors=c("trial", "trialtype"), number.trialtypes=3, iti.mean=0)
print(densityplot(Test_Cued_FPS_mids_Trials, ~FPS)) 
bwplot(Test_Cued_FPS_mids_Trials) # Note, some high values in the original data, less in the imputed data.
# MS & RG: OK

# Test Context FPS --------------------------------------------------------
Test_Context_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag, output_data = out_Umag_trials, M=100, var_pattern = 'sG_Ctx', 
  within.factors=c("trial", "trialtype"), number.trialtypes=3, iti.mean=0)
print(densityplot(Test_Context_FPS_mids_Trials, ~FPS))
bwplot(Test_Context_FPS_mids_Trials)
# MS & RG: OK

# Acquisition ITI (Trials)---------------------------------------------------------
Acquisition_ITI_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag,output_data = out_Umag_trials, M=100, var_pattern = 'sA_ITI_', 
  within.factors=c("trial"), number.trialtypes=0, iti.mean=0) 
print(densityplot(Acquisition_ITI_FPS_mids_Trials, ~FPS ))
bwplot(Acquisition_ITI_FPS_mids_Trials)
# MS & RG: lot of variation between imputations (= quality insufficient), analyze imputed means of ITI (that have good quality imputations).

# Acquisition ITI (Mean)
Acquisition_ITI_FPS_mids_MeanAll <- MAKE_MIDS_Imputed_FPS(
input_data=Umag, output_data=out_Umag_all, M=100, var_pattern = "sA_ITI",
within.factors=c(0), number.trialtypes=0, iti.mean=1)
print(densityplot(Acquisition_ITI_FPS_mids_MeanAll, ~FPS ))
bwplot(Acquisition_ITI_FPS_mids_MeanAll)
# MS & RG: OK

# Test ITI ----------------------------------------------------------------
Test_ITI_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag,output_data = out_Umag_trials, M=100, var_pattern = 'sG_ITI_', 
  within.factors=c("trial"),number.trialtypes=0, iti.mean=0)
# No imputations for this trialtype
# print(densityplot(Test_ITI_FPS_mids_Trials, ~FPS ))
# bwplot(Test_ITI_FPS_mids_Trials)
summary(Test_ITI_FPS_mids_Trials) # for 3 participants, all ITI's were missing. This is more than 2/3, so no imputations were performed.
summary(Umag)

# Test ITI (Mean)
Test_ITI_FPS_mids_MeanAll <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag, output_data=out_Umag_all, M=100, var_pattern = "sG_ITI",
  within.factors=c(0), number.trialtypes=0, iti.mean=1)
print(densityplot(Test_ITI_FPS_mids_MeanAll, ~FPS ))
bwplot(Test_ITI_FPS_mids_MeanAll)
# MS & RG: OK

# pre-Acquisition Habituation ---------------------------------------------
Acquisition_HAB_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag, output_data = out_Umag_trials, M=100, var_pattern = 'sA_HAB_', 
  within.factors=c("trial"),number.trialtypes=0, iti.mean=0)
print(densityplot(Acquisition_HAB_FPS_mids_Trials, ~FPS ))
bwplot(Acquisition_HAB_FPS_mids_Trials)
# MS & RG: OK

# pre-Test Habituation ----------------------------------------------------
Test_HAB_FPS_mids_Trials <- MAKE_MIDS_Imputed_FPS(
  input_data=Umag, output_data = out_Umag_trials, M=100, var_pattern = 'sG_HAB_', 
  within.factors=c("trial"),number.trialtypes=0, iti.mean=0)
print(densityplot(Test_HAB_FPS_mids_Trials, ~FPS ))
bwplot(Test_HAB_FPS_mids_Trials)
# MS & RG: OK


# CREATE OUTPUT -----------------------------------------------------------

# Put elements that can be analyses (so not imputed ITI trials) in a list and save
Umag.mids <-list(
  Acquisition_HAB_FPS_mids_Trials, 
  HABctx_FPS_mids_Trials,
  Acquisition_Cued_FPS_mids_Trials,
  Acquisition_Context_FPS_mids_Trials, 
  #  Acquisition_ITI_FPS_mids_Trials,
  Acquisition_ITI_FPS_mids_MeanAll,
  Test_ITI_FPS_mids_MeanAll,
  #  Test_ITI_FPS_mids_Trials,
  Test_HAB_FPS_mids_Trials,
  Test_Cued_FPS_mids_Trials, 
  Test_Context_FPS_mids_Trials)

names(Umag.mids) <- c(
  "Acq_HAB_Trials", 
  "Acq_HABCTX_Trials", 
  "Acq_Cued_Trials",
  "Acq_Context_Trials", 
  #  "Acq_ITI_Trials",
  "Acq_ITI_Means",
  "Test_ITI_Means",
  #  "Test_ITI_Trials",
  "Test_HAB_Trials", 
  "Test_Cued_Trials",
  "Test_Context_Trials")

save(Umag.mids, file = 'processed_data/Umag.mids.28.01.19.rda') # Sorted Mids objects
