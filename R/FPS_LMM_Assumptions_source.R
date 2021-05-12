#' ---	
#' title: "Check LMM assumption on FPS data from the Fear Generalization Task (FGT) in the SAM study"
#' author: "Milou Sep"	
#' date: '`r format(Sys.time(), "%d %B, %Y")`'
#' output:	
#'   html_document: default	
#' ---	

#' Load required packages	
library(mice); library(lme4); library(influence.ME); library(miceadds)

#' Define required variables
dependent_var='FPS'# Fear Potentiated Startle
n_subjects = 117 # Note this is always the same for all outcome measures (or datasets) in the FGT
influential.case.names<-c()
influentialpoints<-list()

#' Loop assumption checks over all datasets (= list elements) in Data_list
for(i in 1:101){ # Note 101 = 1 dataset with missing values and 100 imputed datasets

# Print in output which dataset is checked --------------------------------
  if (i == 1){print(paste(dataset_name, dependent_var, "- LMER assumption check on data with missings"))
  }else if (i>=1){ print(paste(dataset_name, dependent_var, "- LMER assumption check on imputed dataset", i-1))}  # i-1 is included because the data list contains 1 complete case dataset, so (i=2)-1 is imputed dataset 1

# Select Dataset for assumption check -------------------------------------
  check.data=Data_list[[i]]

# Log-transform FPS (if needed) ----------------------------------------------
  # Note log-transformation needs to be performed prior to the formula call. The unstandardized data contains 0-responses, which will lead to infinite values (that LMER() can not handle). 
  # The solution is to add 1 to all FPS values BEFORE the log-transformation (log()) (ref https://stackoverflow.com/questions/8415778/lm-na-nan-inf-error). 
  # Note, an alternative  solution, to make 0-responses NA, is not suitable here because 0-responses contain meaningful information for FPS analyses.
  # The 0-responses stay in the analyses with +1, as log(1)=0.
  if (log_transformation == 'yes'){
    check.data$FPS <- check.data$FPS + 1
    check.data$FPS <- log(check.data$FPS)
  }
  # print(str(check.data$FPS)) # To check if transformation works.

# Fit Linear(mixed)model --------------------------------------------------
  if (within_factor == 1){ 
    FullModel<-lmer(Model_Formula, data=check.data, REML=F, na.action = na.omit)
  } else if (within_factor == 0){
    FullModel<-lm(Model_Formula, data=check.data)  
  }
  
# Plot's for LMER assumption checks ---------------------------------------
  ## Raw data Histogram
  datacolum<-check.data[which(colnames(check.data) == dependent_var)]
  hist(data.matrix(datacolum), main=paste(dataset_name, "dataset", i, "- Historgram raw data", dependent_var))
  ## Linearity
  plot(fitted(FullModel),residuals(FullModel), main=paste(dataset_name, "dataset", i, "- Linearity contrast", dependent_var)); abline(0,0,lty=2,lwd=2);
  ## Normality of residuals
  ### histogram Residuals
  hist(residuals(FullModel), main=paste(dataset_name, "dataset", i, "- Historgram Model Residuals", dependent_var));
  ## QQ-plots
  qqnorm(residuals(FullModel), main=paste(dataset_name, "dataset", i, "- QQ plot", dependent_var)); qqline(residuals(FullModel),col = 2)
  ## Cook's distance
  if (within_factor == 1){ # Display Cook's distance for LMER() model
    estex.FullModel <-influence(FullModel,"pp")
    cooks_FullModel<-cooks.distance(estex.FullModel)
    these_influential_points<-cooks_FullModel[(cooks_FullModel>(4/(n_subjects))),] # Shows measures with cook's distance above cutoff
    print(these_influential_points)
    # Plot cook's distance
    plot(estex.FullModel, which="cook",
         cutoff=4/(n_subjects), sort=TRUE,
         xlab="Cook's Distance",
         ylab="pp",
         cex.lab=0.01, # Smaller font size for as labels.
         main=paste(dataset_name, "dataset", i, "Cook's Distance", dependent_var))
    # Collect influential cases (as information for analyses)
    influentialpoints[[i]] <- these_influential_points  # participant name & value of cook's distance
    if (i>1){ # collects only names of participants that are influential in imputed datasets(not cooks distance data)
      influential.case.names <- c(influential.case.names, names(influentialpoints[[i]]))
    }
  } else if (within_factor == 0){ # Display Cook's distance for lm() model
    plot(FullModel,2) # QQ plot, with case names.
    plot(FullModel,4) # Plot for cooks distance for lm() (Note, this shows observations, not participants numbers)
  }
}

#' Remove redundant variables from the global environment.
rm(list=c("check.data", "cooks_FullModel", "Data_list", "datacolum", "dataset_name", "dependent_var", "estex.FullModel", 
          "FullModel", "i", "Model_Formula", "n_subjects", "these_influential_points", "within_factor"), pos=.GlobalEnv) 
