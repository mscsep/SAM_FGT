## Script to calls "FPS_LMM_Assumptions_scours.r" multiple times, to check the assumptions for the linear(mixed)model analyses in the data from the Fear Generalization Task (FGT) in the SAM study
# Written by Milou Sep

## Selection of Random effects in full models:
# Although, factors 'Trial' and 'Trialtype' are both 'within factors' in most datasets below. Random slopes are only added to the full model for subject (+1|pp).
# To compare nested models, random effects need to be the same, therefore all (LMER) models have random effects: +(1|pp). Various random effects were checked, 
# but more complex random effects than 1|pp (like 1|pp+1|trialtype or 1|pp+1|trial) cause model convergence problems in some outcome measures. 
# We therefore selected the simplest random effects term (note various random effects did not crucially influence model results, i.e. no difference in the effects that reached significance.)

# Set WD & Load Mids ------------------------------------------------------
load("processed_data/Umag.mids.28.01.19.rda")
summary(Umag.mids)

# Description of the variables that are defined prior to each 'assumption check'. Assumptions are checked by rendering the file `FPS_LMM_Assumptions_source.R`.
# dataset_name --> just a character string that is used in the .html output
# Data_list --> a list of datasets. The 1 element is the dataframe with missing values, the following are the 100 imputed datasets
# log_transformation --> can be "yes" or something else. If "yes", a log-transformation of FPS is performed in the source script
# Model_Formula=formula(" ") --> the formula for LM() or LMER() should be added as a character string
# within_factor --> can be 1 (=present) or 0 (=absent). Cooks distances are differently calculated for LM() (if no within factor is present) and lmer() (if one or multiple within factors are present)

###### ACQUISITION ######

# Habituation Trials Day 1 -------------------------------------------------
# Required variables definitions
{dataset_name = 'Hab_trials_D1'
Data_list=c(list(Umag.mids$Acq_HAB_Trials$data),
            miceadds::mids2datlist(Umag.mids$Acq_HAB_Trials, X=NULL))
log_transformation = 'yes'
Model_Formula=formula("FPS~1+Condition*trial+(1|pp)") 
within_factor = 1
# Render Assumption Check script to HTML 
rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document", 
                  output_file = paste0('Assumptions_HAB_trials_D1.',date(),'.html')) 
}# Show individual influential points
save(influential.case.names, file=paste0('influentialCases_HAB_trials_D1.',date(),'.rda'))
influential.case.names # Show influential participants:
unique(influential.case.names)
mean(unlist(influentialpoints), na.rm=T) # Mean cook's distance of outliers.
rm("influentialpoints", "influential.case.names", "log_transformation")


# Cue Trials Day 1 -------------------------------------------------
# Required variable definitions
{dataset_name = 'Cue_trials_D1'
Data_list=c(list(Umag.mids$Acq_Cued_Trials$data),
            miceadds::mids2datlist(Umag.mids$Acq_Cued_Trials, X=NULL))
log_transformation = 'yes'
Model_Formula=formula("FPS~1+Condition*trial*trialtype+(1|pp)") 
within_factor = 1
# Render Assumption Check script to HTML 
rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document", 
                  output_file = paste0('Assumptions_Cue_trials_D1.',date(),'.html'))
}# Show individual influential points
save(influential.case.names, file=paste0('influentialCases_Cue_trials_D1.',date(),'.rda'))
influential.case.names 
#unique(influential.case.names)
mean(unlist(influentialpoints), na.rm=T) # Mean cook's distance of outliers.
rm("influentialpoints", "influential.case.names", "log_transformation")


# Context Trials Day 1 -------------------------------------------------
# Required variable definitions
{dataset_name = 'Context_trials_D1'
Data_list=c(list(Umag.mids$Acq_Context_Trials$data),
            miceadds::mids2datlist(Umag.mids$Acq_Context_Trials, X=NULL))
log_transformation = 'yes'
Model_Formula=formula("FPS~1+Condition*trial*trialtype+(1|pp)") 
within_factor = 1
# Render Assumption Check script to HTML 
rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document", 
                  output_file = paste0('Assumptions_Context_trials_D1.',date(),'.html'))
}# Show individual influential points
save(influential.case.names, file= paste0('influentialCases_Context_trials_D1.', date(),'.rda'))
influential.case.names 
unique(influential.case.names) 
mean(unlist(influentialpoints), na.rm=T) # Mean cook's distance of outliers.
rm(influentialpoints, influential.case.names, "log_transformation")


# ITI Mean Day 1 -------------------------------------------------
# Required variable definitions
{dataset_name = 'ITI_Mean_D1'
Data_list=c(list(Umag.mids$Acq_ITI_Means$data),
            miceadds::mids2datlist(Umag.mids$Acq_ITI_Means, X=NULL))
log_transformation = 'yes'
Model_Formula=formula("FPS~1+Condition") 
within_factor = 0
# Render Assumption Check script to HTML 
rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document", 
                  output_file = paste0('Assumptions_ITI_Mean_D1.',date(),'.html'))
} # Note variables "influential.case.names" and "influential points" are only created when LMER() is used (in this case LM() is used, because there is only 1 grouping factor).
rm(influentialpoints, influential.case.names, "log_transformation")



###### TEST ######

# Habituation Trials Day 2 -------------------------------------------------
# Required variable definitions
{dataset_name = 'Hab_trials_D2'
  Data_list=c(list(Umag.mids$Test_HAB_Trials$data),
              miceadds::mids2datlist(Umag.mids$Test_HAB_Trials, X=NULL))
  log_transformation = 'yes'
  Model_Formula=formula("FPS~1+Condition*trial+(1|pp)") 
  within_factor = 1 
  # Render Assumption Check script to HTML 
  rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document", 
                    output_file = paste0('Assumptions_HAB_trials_D2.',date(),'.html'))
}# Show individual influential points
save(influential.case.names, file= paste0('influentialCases_HAB_trials_D2.',date(),'.rda'))
influential.case.names 
unique(influential.case.names)
mean(unlist(influentialpoints), na.rm=T) # Mean cook's distance of outliers. 
rm("influentialpoints", "influential.case.names", "log_transformation")


# Cue Trials Day 2 -------------------------------------------------
# Required variable definitions
{dataset_name = 'Cue_trials_D2'
Data_list=c(list(Umag.mids$Test_Cued_Trials$data),
            miceadds::mids2datlist(Umag.mids$Test_Cued_Trials, X=NULL))
log_transformation = 'yes'
Model_Formula=formula("FPS~1+Condition*trial*trialtype+(1|pp)")
within_factor = 1
# Render Assumption Check script to HTML
rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document",
                  output_file = paste0('Assumptions_Cue_trials_D2.',date(),'.html'))
}# Show individual influential points
save(influential.case.names, file= paste0('influentialCases_Cue_trials_D2.',date(),'.rda'))
influential.case.names
unique(influential.case.names) 
length(influential.case.names)
mean(unlist(influentialpoints), na.rm=T) # Mean cook's distance of outliers.
rm(influentialpoints, influential.case.names, "log_transformation")


# Context Trials Day 2 -------------------------------------------------
# Required variable definitions
{dataset_name = 'Context_trials_D2'
Data_list=c(list(Umag.mids$Test_Context_Trials$data),
            miceadds::mids2datlist(Umag.mids$Test_Context_Trials, X=NULL))
log_transformation = 'yes'
Model_Formula=formula("FPS~1+Condition*trialtype*trial+(1|pp)")
within_factor = 1
# Render Assumption Check script to HTML
rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document",
                  output_file = paste0('Assumptions_Context_trials_D2.',date(),'.html'))
}# Show individual influential points
save(influential.case.names, file= paste0('influentialCases_Context_trials_D2.', date(), '.rda'))
influential.case.names
unique(influential.case.names) 
length(influential.case.names)
mean(unlist(influentialpoints), na.rm=T) # Mean cook's distance of outliers.
rm(influentialpoints, influential.case.names, "log_transformation")


# ITI Mean Day 2 -------------------------------------------------
# Required variable definitions
{dataset_name = 'ITI_Mean_D2'
  Data_list=c(list(Umag.mids$Test_ITI_Means$data),
              miceadds::mids2datlist(Umag.mids$Test_ITI_Means, X=NULL))
  log_transformation = 'yes' # Log-transformation is needed, to ensure normal distribution of the residuals.
  Model_Formula=formula("FPS~1+Condition") 
  within_factor = 0
  # Render Assumption Check script to HTML 
  rmarkdown::render("R/FPS_LMM_Assumptions_source.R", "html_document", 
                    output_file = paste0('Assumptions_ITI_Mean_D2.',date(),'.html'))
} # Note variables "influential.case.names" and "influential points" are only created when LMER() is used (in this case LM() is used, because there is only 1 grouping factor).
rm(influentialpoints, influential.case.names, "log_transformation")
