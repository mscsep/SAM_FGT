# Imputation Fear potentiated startle (EMG eyeblink) data from the Fear Generalization Task (FGT) in the SAM study
# Written by Rosalie Gorter & Milou Sep.

# NOTE: The function "Impute_FGT_EMG_SAM" is used twice in this script: 
# 1) to impute individual trials (set sorttype to 'trials'). 
#     - Note: trials will be imputed if more than 1/3 of the trials in a category is present (in other words if <2/3 missing), 
#       if less than 1/3 is present (in other words if >2/3 is missing; `missing code 4`) all the trials (for that category) will be set to missing 
# 2) to create imputed means (set sorttype to 'mean')
#     - Note: means will be based on imputed trials if more than 2/3 of trials is present (in other words if <1/3 missing; `missing code 1`), 
#      or imputed directly if less than 2/3 of the trials is present (or in other words if >1/3 missing; `missing code 2` ).

#install.packages("mice")
library("mice")
# install.packages("readr")
library(readr)

# Read Data ---------------------------------------------------------------
# (Note: Only valid digital Biopac cues analyzed)
FGT_batch  <- read_delim("data/SAM_FGT.csv", ";",locale=locale(decimal_mark = ","), escape_double = FALSE, trim_ws = TRUE, na = c("NaN","5555","8888","9999"))
# Note Missing codes (assigned via Matlab preprocessing code):
# 5555 % Digital registration error, message on screen & MissingValue 5555.
# 8888 % Excessive baseline activity. (More than 2 SD deviation)
# 9999 % Latency onset not valid (Note, this only affects onset scores, so not present in magnitude variables)

# Prepare Data ------------------------------------------------------------
# remove "SMGH" from subject names to use numbers on x-as (& in subset)
FGT_batch$subjName <- gsub("SMGH","", FGT_batch$subjName)
# Make Condition Factor
FGT_batch$Condition <- as.factor(FGT_batch$Condition)
# Subset all Unstandardized startle magnitude variables (Umag)
Umag <- subset.data.frame(FGT_batch, select = c(grep("Umag", names(FGT_batch)),Condition))
# Note! Imputation now only done on Umag data (which is better than standardized data for LMM analyses). The dataset also contains Standardized startle magnitude variables (Smag), Peak latencies () and Onset Latencies are also available in FGT_batch. 
# This script could also be used to impute Standardized FPS responses.
# Smag <- subset.data.frame(FGT_batch, select = c(grep("Smag", names(FGT_batch)),Condition))

# Imputation Settings -----------------------------------------------------
M<-100
MAXIT<-50

# Call imputation script --------------------------------------------------

# Impute separate trials (used for Habituation-trials, Cue-trials, Context-trials)
out_Umag_trials<-Impute_FGT_EMG_SAM(data=Umag, M=M, MAXIT=MAXIT, sorttype='trials') # call function below
save(out_Umag_trials, file = 'processed_data/OUPUT_IMPUTATIE_FGT_out_M50_MAXIT100_Umag_Trials.rda') # save Imputed datasets

# Impute means (used for inter-trial intervals)
out_Umag_all<-Impute_FGT_EMG_SAM(data=Umag, M=M, MAXIT=MAXIT, sorttype='mean') # call function below
save(out_Umag_all, file = 'processed_data/OUPUT_IMPUTATIE_FGT_out_M50_MAXIT100_Umag_AllMeans.rda') # save Imputed datasets

# Export original data ------------------------------------------------------
# Add subjectname to dataset (was removed for imputations, but required for further data analyses)
Umag$subject <- as.numeric(FGT_batch$subjName)
# Save input data to imputation file
save(Umag, file = 'processed_data/INPUT_IMPUTATIE_FGT_Umag.rda') # Raw dataset (needed to select `complete cases` for midsobject)


# Imputation Function Definition ------------------------------------------
Impute_FGT_EMG_SAM<-function(data,M,MAXIT,sorttype,...){
  
  # DEFINE VARIABLE PATTERNS & NAMES ----------------------------------------
  FGT_Variabels <- c('sA_Cue_.*_T', 'sA_Cue_.*_S', 'sA_Ctx_.*_T','sA_Ctx_.*_S', 'sA_ITI_',   # Variables 1t/m5
                     'sA_HAB_', # Variable 6
                     'sA_HABctx_.*_T', 'sA_HABctx_.*_S','sA_HABctx_.*_N', # Variables 7,8,9
                     'sG_Cue_.*_T', 'sG_Cue_.*_S', 'sG_Cue_.*_N',   # Variables 10,11,12
                     # Variables with only 3 trials per category:
                     'sG_Ctx_.*_T','sG_Ctx_.*_S', 'sG_Ctx_.*_N', 'sG_ITI_', # Variable 13, 14, 15, 16
                     'sG_HAB_' ) # Variable 17 [same method as variable 6]
  FGT_Variabels_Names <- gsub(x=FGT_Variabels, pattern="_.*_", replacement='_', fixed=TRUE) # Note Fixed true is needed for r to 'interpret' .* as character and not as 'wildcard'
  
  # SELECT & RESTRUCTURE FOR IMPUTATION -------------------------------------
  restructure_FGTdata <-function(data_restruct, sorttype){ # Function to restructure data
    # Make empty lists for the for-loop below
    selected_FGT_data_allTrials <- list()
    # In de for-loop below are all columns selected
    # Select the columns that fit the variables in FGT_Variables vector
    for(i in 1:length(FGT_Variabels)){
      this_data <- data_restruct[,(grep(FGT_Variabels[i], names(data_restruct)))]
      selected_FGT_data_allTrials[[i]] <- this_data
    }
    # Name elements in lists
    names(selected_FGT_data_allTrials) <- FGT_Variabels_Names
    # Return variable
    return(selected_FGT_data_allTrials)
  }
  
  ## Call restructure script that is defined above ##
  ldata <- restructure_FGTdata(data, sorttype)
  
  
  # COUNT & RECODE MISSINGS -------------------------------------------------
  Recode_Missings_FGT<-function(sorttype){
    # Function to count & code missing values in FGT startle data. 
    
    # 1) Specific for the imputation of means (sorttype = "all"):
    # Note: Codes 0 = All present (0 missing values), 1 = Less than 1/3 missing ((or in other words if more than 2/3 present), 2 = More than 1/3 missing (or in other words if less than 2/3 present) 
    # Result: Code 1 creates means based on imputed trials, Code 2 indicates that the mean needs to be imputed directly.
    if (sorttype == 'mean') { # If more than 1/3 missing: Code 2 
      n.missings <- code.missings <- list() # make empty list
      for (i in 1:length(ldata)) {
        n.missings[[i]] <- rowSums(is.na(cbind(ldata[[i]][, 1:ncol(ldata[[i]])])  * 1)) # Note n.missing is the number of missing trials PER participant WITHIN one category (e.g. Acquisition threat)
        code.missings[[i]] <-
          ifelse(n.missings[[i]] > (1/3 * ncol(ldata[[i]])), 2, 1) # [Note: ifelse(test, yes, no)]
        code.missings[[i]] <-
          ifelse(n.missings[[i]] == 0, 0, code.missings[[i]])}
      
      # 2) Specific for the imputation of individual trials (sorttype = 'trials'):
      # Note: Code 0 = All present (0 missing values); Code 4 = More than 2/3 missing (or in other words if less than 1/3 present)
      # Result: Code 4 indicates that all trials will be set to missing; without code 4 (so less than 2/3 missing, or in other words more than 1/3 present) the individual trials are imputed.
    } else if (sorttype == 'trials'){ 
      n.missings <- code.missings <- list()
      for (i in 1:length(ldata)) {
        n.missings[[i]] <- rowSums(is.na(cbind(ldata[[i]][, 1:ncol(ldata[[i]])])  * 1))  
        code.missings[[i]] <-
          ifelse(n.missings[[i]] > (2/3 * ncol(ldata[[i]])), 4, 1) # If more than 2/3 missing: Code 4
        code.missings[[i]] <-
          ifelse(n.missings[[i]] == 0, 0, code.missings[[i]])}
    }
    
    out <- list(code.missings, n.missings) 
    names(out) <- c("code", "n")
    return(out)
  }
  
  # Call function above. Note: `Missings` is a list of two with `Missings[[1]]` = code.missings and `Missings[[2]]` = n.missings
  Missings <- Recode_Missings_FGT(sorttype)
  
  # Save trials & missings in list and as file
  data_missings_trials<-list(ldata,Missings$n,Missings$code)
  names(data_missings_trials)<-c("ldata","n.missings","code.missings")
  saveRDS(data_missings_trials, "processed_data/data_missing_trials.rda")
  
  
  if(sorttype == 'mean'){ # If imputation should create imputed means:
    
    # CALCULATE MEAN STARTLES -------------------------------------------------
    mean.startle<-matrix(unlist(lapply(ldata,rowMeans,na.rm=T)),nrow=nrow(data)) # True = omit missing values van calculations
    mean.startle.completecases<-matrix(unlist(lapply(ldata,rowMeans,na.rm=F)),nrow=nrow(data)) # False = do not omit missing values van calculations [result: if a missing value is present, the mean is missing]
    # Create & Assign variable names for means
    FGT_Variabels_Names_Means <- paste0("mean_",FGT_Variabels_Names)
    colnames(mean.startle) <-  FGT_Variabels_Names_Means
    colnames(mean.startle.completecases) <-  FGT_Variabels_Names_Means
    
    # SAVE MEANS --------------------------------------------------------------
    #save data in list and in file
    data_missings_means<-list(ldata,Missings$n,Missings$code,mean.startle)
    names(data_missings_means)<-c("ldata","n.missings","code.missings","mean.startle")
    saveRDS(data_missings_means, "processed_data/data_missings_means_all.rda")
    
    # PREPARE DATA FOR IMPUTATION ---------------------------------------------
    # generate dataset with only means based on more than 2/3 of trials (pulses)
    for(i in 1:length(ldata)){
      for(j in 1:nrow(data)){
        if(Missings$code[[i]][j]==2){mean.startle[j,i]<-NA} # change the means that are based on less than 2/3 of trials (pulses) (code 2) to NA
      }
    }
    # Add Condition to the matrix with means
    mean.startle.c<-cbind(mean.startle,data$Condition)
    colnames(mean.startle.c)[length(colnames(mean.startle.c))]<-"Condition"
    
  }
  
  
  # MAKE PREDICTION MATRIX --------------------------------------------------
  # Note: the imputation on all trials is performed before trials are sorted in categories (so applicable to sorttype = 'trials' and sorttype = 'mean')
  pred_allpulses = (1 - diag(1, ncol(data)))
  rownames(pred_allpulses) <- colnames(data)
  colnames(pred_allpulses) <- colnames(data)
  pred_allpulses["Condition", 1:ncol(pred_allpulses)] <- 0
  pred_allpulses[1:ncol(pred_allpulses)-1, "Condition"] <- 1
  
  # Prediction matrix for imputation of means: specific for sorttype = 'mean' (not needed if sorttype = 'trials').
  if (sorttype == 'mean'){
    pred_means = (1 - diag(1, ncol(mean.startle.c)))
    pred_means[18, 1:ncol(pred_means)] <- 0 # The number (18) indicate the variable location + 1 (because condition was added to matrix)
    pred_means[1:ncol(pred_means)-1, 18] <- 1
  }
  
  # IMPUTATION --------------------------------------------------------------
  ## Trials
  Startle_imputatie_pulses <- mice(data=data, 
                                   pred=pred_allpulses, 
                                   m = M, 
                                   maxit = MAXIT, 
                                   seed=356)
  saveRDS(Startle_imputatie_pulses, "startle_pulses_imp")
  ## Means
  if (sorttype == 'mean'){
    Startle_imputatie_means <- mice(mean.startle.c, 
                                    pred=pred_means, 
                                    m = M, 
                                    maxit = MAXIT, 
                                    seed=356)
    saveRDS(Startle_imputatie_means, "startle_means_imp")
  }
  
  # MERGE DATA --------------------------------------------------------------
  dataout <- list() # Empty list
  for(m in 1:M) { # number of Imputed datasets
    
    ## Create 3 complete datasets, and add suffix _i to all variable names:
    
    # 1) Dataset with imputed trials (pulses) (for mean calculations & individual trials dataset)
    Imputed_Pulses_complete <- complete(Startle_imputatie_pulses, m,include = F) # Note include = T, would include original data with missing values included.
    colnames(Imputed_Pulses_complete) <- paste0(names(Imputed_Pulses_complete),"_i")
    
    # 2) Only when sorttype = 'mean'
    if (sorttype == 'mean'){
      # 2A) Dataset with imputed means
      Imputed_Means_complete <- complete(Startle_imputatie_means, m,include = F)
      colnames(Imputed_Means_complete) <- paste0(names(Imputed_Means_complete),"_i")
      # 2B) Dataset with calculated means based on dataset 1: the imputed trials (pulses)
      ldataIMP <- restructure_FGTdata(Imputed_Pulses_complete, sorttype) # Imputed data sorted according to sorttype, using the 'restructure function' 
      Means_based_on_Imputed_Pulses<-list() # Create empty list
      for (i in 1:length(ldataIMP)) { # fill empty mean list
        this_mean <- rowMeans( cbind(ldataIMP[[i]]) ,na.rm=TRUE) # Calculate row means for all elements in ldataIMP
        Means_based_on_Imputed_Pulses <- as.data.frame(cbind(Means_based_on_Imputed_Pulses, this_mean))} # add row mean to dataframe
      # Add column names to dataset 2B
      if (sorttype == 'mean'){colnames(Means_based_on_Imputed_Pulses) <- paste0(FGT_Variabels_Names_Means,"_i")}
      ## Create Merged means dataset. Note this dataset contains the appropriate mean (either from dataset 2A or 2B), based on the Missing$Code.
      mergedmeans <- data.frame (matrix(0,nrow(data), length(ldata))) # Make empty matrix
      # Add column names + suffix 'merged'
      if (sorttype == 'mean'){colnames(mergedmeans) <- paste0(FGT_Variabels_Names_Means,"_merged")}
      # Fill the `mergedmeans` matrix with the appropriate data, based on missing codes.
      for(i in 1:length(ldata)){
        for(j in 1:nrow(data)){
          if(Missings$code[[i]][j]==1){mergedmeans[j,i]<-Means_based_on_Imputed_Pulses[j,i]} # if less than 1/3 Missing (or 1/3 missing)
          if(Missings$code[[i]][j]==0){mergedmeans[j,i]<-Imputed_Means_complete[j,i]} # 0 missing
          if(Missings$code[[i]][j]==2){mergedmeans[j,i]<-Imputed_Means_complete[j,i]} # if more than 1/3 missing
        }
      }
      dataout[[m]] <- mergedmeans
      
      # 3) Only when sorttype = 'trials'
    } else if (sorttype == 'trials'){ 
      # 3A) dataset with separate imputed trials, corrected for the `2/3 missing rule`, for analyses with "trials" as factor.
      # Restructure imputed pulses (trials) (of imputed dataset m)
      ldataIMP <- restructure_FGTdata(Imputed_Pulses_complete, sorttype) # Imputed data sorted according to sorttype, using the 'restructure function' (output in list)
      # copy data to new list before transformations
      ldataIMP_corrected <- ldataIMP
      # Loop over all variable categories (the different elements in ldataIMP, 17 in total), e.g. threat acquisition, safe acquisition etc.
      for (i in 1:length(ldataIMP)){   # length ldataIMP = 17 
        # Loop per variable category over all participants 
        for (j in 1:nrow(data)){ # nrow(data) = 117
          # Check per variable category i for each participant j if the data (in that variable category) needs recoding.
          # Missing Codes: 0 = All present, 4 = More than 2/3 missing
          if(Missings$code[[i]][j]==4){ # when more than 2/3 of the trials in this category is missing
            ldataIMP_corrected[[i]][j,] <-  rep(NA,length(ldataIMP[[i]])) # rownumber j in dataset i, replaced by NA's for all trials (= length ldataIMP[[i]])
            # For detailed checking of missing codes (if required)
            # Missings$code[[1]]# a vector with a missing code for each participant in Threat acquisition (dataset 1) trials
            # Missings$code[[1]][4] # a vector with a missing code for participant 4 in Threat acquisition (dataset 1) trials
          }
        } 
      }
      
      # 3B Convert list, to dataframe
      All.Trials.MissingCorrected<-cbind(data.frame(matrix(unlist(ldataIMP_corrected),nrow = 117)), Imputed_Pulses_complete$Condition_i)
      colnames(All.Trials.MissingCorrected) <- colnames(Imputed_Pulses_complete)  # add names
      # Store in output list
      dataout[[m]]<-All.Trials.MissingCorrected
    }
    
  }
  
  
  # CREATE OUTPUT -----------------------------------------------------------
  if (sorttype != 'trials'){
    out<-list(dataout,mean.startle.completecases)
    names(out)<-c("MergedMeans_imputed","Means_Completecases")
  } else if (sorttype == 'trials'){
    out<-list(dataout, data)
    names(out)<-c("Trials_imputed","Trials_Original")
  }
  
  # Save output in the processed_data folder
  save(out, file = 'processed_data/OUPUT_IMPUTATIE_FGT_out_endofscript.rda') # Imputed datasets
  
  return(out)
  
}
