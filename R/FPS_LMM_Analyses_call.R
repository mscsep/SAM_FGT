#' ---	
#' title: "LMM Analyses of unstandardized Fear Potentiated Startle responses to Cue, Context, ITI and Habituation probes (Learning D1 and Memory D2) in the Fear Generalization Task (FGT) in the SAM study"
#' author: "Milou Sep"	
#' output:	
#'   html_document: default	
#' ---	

# some information on analyses of imputed data: http://web.maths.unsw.edu.au/~dwarton/missingDataLab.html
#' # FGT "fear conditioning" specific analysis info: 
#+ Section, include = FALSE
#   * separate tasks phases (separate analyses)
#     - day 1 (learning) and day 2 (memory)
#   * outcome measures in FGT (separate analyses):
#     - startle probe habituation (=hab)
#     - context habituate [only day 1] (=habctx, Trialtype: threat, safe & new)
#     - cue & context learning/memory (=cue and context; Trialtype: threat, safe & new; Startle probetypes: cue vs context)
#     - inter-trial interval (=ITI)

# Loading required packages -----------------------------------------------
library(mice); library(lme4); library(miceadds)

# Loading imputed data ----------------------------------------------------
load("processed_data/Umag.mids.28.01.19.rda") # imputed data: 1) combined Mids object, and day 1 & day 2 Mids subsets. 
# Note These files contain separate sorted Mids objects for cue and context, habctx, hab (per day). Note habctx is only measured on day1.	
# Note For ITI on day 1 and day 2 are the imputed means analyzed (because imputation quality of individual trials was insufficient, quality of the imputed means was good).

# Source analyses functions -----------------------------------------------
source("R/FPS_LMM_trialtype.trial.condition.r")
source("R/FPS_LMM_trialtype.condition.r")
source("R/FPS_LMM_trial.condition.R") # For analyzes of responses to Habituation probes (on d1 and d2)
source("R/FPS_LM_condition.r") # For analyses mean ITI (D1 & D2)
source("R/FPS_LMM_pool_EMM.r") # Contains 1) RG's function to pool LMER estimated marginal means (EMM) from imputed datasets (used to followup "Trialtype" Effects in early-mid-late epochs) and 2) a function to plot mean epochs of follow-up analyses
# Load data export function
source("R/FPS_LMM_Results_to_Word.R") # Function to export results in word-friendly format
# Source functions for data transformations mids objects
source("R/FGT_mids_transformations.R")

# Sensitivity analyses with participants that detected the threat stimulus correctly (according to self-report). Note this information was not collected for the complete sample.
# These sensitivity analyses were just for exploration (not published).
readRDS("processed_data/participants.for.contingency.sensitivity.analyses.rda")->pp_Sensitivity # Data in processed_data folder (created with "FGT_descriptive.R")
pp_Sensitivity$SubjectNumber <- gsub("SMGH","", pp_Sensitivity$SubjectNumber) # Remove "SMGH" before participant number
pp_Sensitivity$SubjectNumber <-as.factor(as.numeric(pp_Sensitivity$SubjectNumber)) # Change to numeric and than to factor (=same as mids)
str(pp_Sensitivity$SubjectNumber) # to check


# Call of Analyses scripts ------------------------------------------------

# Day 1: Learning / Acquisition phase -------------------------------------

# Habituation (Noise Alone) FPS (acquisition) ---------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
A_Hab_Trials.log <- log.transform.mids(Umag.mids$Acq_HAB_Trials)
Acq_Hab_2way <- FGT_LMER_trial.Condition(dataset=A_Hab_Trials.log)
Export.FGT.Table(FGT_Results_Table = Acq_Hab_2way, TableName = "Umag.Hab.Acq.Trials", file.export = T)  # Note this function is in 'FGT_Results_to_Word.R'
# Main effect trial, no interactions (No follow-up analyses planned/needed).

# Sensitivity analyses to check the effect of potential influential participants, identified via Assumption Checks: pp 98 and pp 100
A_Hab_Trials.log.noinf <- remove.influential.points(data_mids=A_Hab_Trials.log, influentialpoints=c(98,100))
Acq_Hab_2way_noinf <- FGT_LMER_trial.Condition(dataset=A_Hab_Trials.log.noinf)
Export.FGT.Table(FGT_Results_Table = Acq_Hab_2way_noinf, TableName = "Umag.Hab.Acq.noinf.Trials", file.export = F)
# Conclusion: no significant change in result.


# Cued FPS (acquisition) --------------------------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
A_Cue_Trials.log <- log.transform.mids(Umag.mids$Acq_Cued_Trials)
Acq_Cued_3way <- FGT_LMER_Trialtype.Trial.Condition(dataset=A_Cue_Trials.log)
Export.FGT.Table(FGT_Results_Table = Acq_Cued_3way, TableName = "Umag.Cue.Acq.Trials", file.export = T)
# Interaction effect Trialtype * Trial (and main effect trialtype & trial)

# Sensitivity analyses to check the effect of potential influential participants, identified via Assumption Checks: pp 15, 97 (and 90?)
A_Cue_Trials.log.noinf <- remove.influential.points(data_mids=A_Cue_Trials.log, influentialpoints=c(15,97))
Acq_Cued_3way_noinf <- FGT_LMER_Trialtype.Trial.Condition(dataset=A_Cue_Trials.log.noinf)
Export.FGT.Table(FGT_Results_Table = Acq_Cued_3way_noinf, TableName = "Umag.Cue.Acq.noinf.Trials", file.export = F)
# Conclusion: no significant change in result.

# * Follow-up analyses: decompose factor "Trial" by analyzing means of early & late trials separately
A_Cue<-complete(A_Cue_Trials.log, action = "long", include = T) # Change mids to long format
unique(A_Cue$trial) # 12 trials

# Create Early Trials mids object
A_Cue_Early <- subset(A_Cue, trial %in% c(1,2,3,4)) 
unique(A_Cue_Early$trial) # to check
A_Cue_Early_M <- means_EML_2way(A_Cue_Early)
# Analyze Early Trials
Acq_Cued_2way.Early <- FGT_LMER_trialtype.Condition(dataset=A_Cue_Early_M)
Export.FGT.Table(FGT_Results_Table = Acq_Cued_2way.Early, TableName = "Umag.Cue.Acq.Early", file.export = T)
# No significant effects. 

# Create Mid Trials mids object
A_Cue_Mid <- subset(A_Cue, trial %in% c(5,6,7,8)) 
unique(A_Cue_Mid$trial)
A_Cue_Mid_M <- means_EML_2way(A_Cue_Mid)
# Analyze Mid Trials
Acq_Cued_2way.Mid <- FGT_LMER_trialtype.Condition(dataset=A_Cue_Mid_M)
Export.FGT.Table(FGT_Results_Table = Acq_Cued_2way.Mid, TableName = "Umag.Cue.Acq.Mid", file.export = T)
# No significant effects. 

# Create Late Trials mids object
A_Cue_Late <- subset(A_Cue, trial %in% c(9,10,11,12)) 
unique(A_Cue_Late$trial)
A_Cue_Late_M <- means_EML_2way(A_Cue_Late)
# Analyze Late Trials
Acq_Cued_2way.Late <- FGT_LMER_trialtype.Condition(dataset=A_Cue_Late_M)
Export.FGT.Table(FGT_Results_Table = Acq_Cued_2way.Late, TableName = "Umag.Cue.Acq.Late", file.export = T)
# Significant effect Trialtype

# * Follow-up Main effect Trialtype (to determine direction of effect): ####
# To  date,  there  are  no statistical methods implemented in the major software packages to pool the results from post-hoc tests on EMMs of imputed datasets. 
# Therefore, differences in EMMs could not be tested statistically in the current study. Also see https://www.tandfonline.com/doi/abs/10.1080/00273171.2013.855890
# The EMM were pooled to get an indication for the direction of the effect.

# Copy Relevant Model from LMER script
A_Cue_Main_NO_2way<-with(expr=lmer(FPS~1+Condition+trialtype+(1|pp), control=lmerControl(optimizer="bobyqa")),data = A_Cue_Late_M, REML=F)
# Call function to pool EMM:
Acq_Cued_2way.Late.followupTrialtype <-pool.means.trialtype(model=A_Cue_Main_NO_2way, m=100, phase="ACQ") # Enter larger model from the specific comparison that is followed up. In this case Full main effects model.
Acq_Cued_2way.Late.followupTrialtype # Note Qbar = pooled means (overall point estimate)
# Export results to table:
Export.FGT.FollowupTrialtype(FGT_FollowupResults = Acq_Cued_2way.Late.followupTrialtype, 
                             TableName = "Umag.Cue.Acq.Late_FU.trialtype", file.export = T)
# No further follow-up analyses planned/needed


# Contextual FPS (acquisition) --------------------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
A_Ctx_Trials.log <- log.transform.mids(Umag.mids$Acq_Context_Trials)
# A_Ctx_Trials.log<-Sensitivity_contingency(A_Ctx_Trials.log) # Same results
Acq_Ctx_3way <- FGT_LMER_Trialtype.Trial.Condition(dataset=A_Ctx_Trials.log)
Export.FGT.Table(FGT_Results_Table = Acq_Ctx_3way, TableName = "Umag.Ctx.Acq.Trials", file.export = T) 
# Main trialtype & trial

# No Sensitivity analyses to check effect influential participants: No influential points were identified

# * Follow-up analyses 1: overall mean effect trialtype (Note these EMM are averaged over Conditions & Trials!)
Main_NO_2ways_ACQCTX<-with(expr=lmer(FPS~1+Condition+trialtype+trial+(1|pp), control=lmerControl(optimizer="bobyqa")),data = A_Ctx_Trials.log, REML=F)
ACQ_CTX_2way.followupTrialtype<- pool.means.trialtype(model=Main_NO_2ways_ACQCTX, m=100, phase='ACQ')
Export.FGT.FollowupTrialtype(FGT_FollowupResults=ACQ_CTX_2way.followupTrialtype, 
                             TableName = "Umag.Ctx.ACQ.ALL_FU.trialtype", file.export = T)
# The effects plotted
Trialtype.PlotD2(ACQ_CTX_2way.followupTrialtype, rownames(ACQ_CTX_2way.followupTrialtype), ACQ_CTX_2way.followupTrialtype$Qbar)

# * Follow-up analyses 2: decompose factor "Trial" by analyzing means of early & late trials separately
A_Ctx<-complete(A_Ctx_Trials.log, action = "long", include = T) # Change mids to long format
unique(A_Ctx$trial) # 6 trials

# Create Early Trials mids object
A_Ctx_Early <- subset(A_Ctx, trial %in% c(1,2)) 
unique(A_Ctx_Early$trial)
A_Ctx_Early_M <- means_EML_2way(A_Ctx_Early)
# Analyze Early Trials
Acq_Ctx_2way.Early <- FGT_LMER_trialtype.Condition(dataset=A_Ctx_Early_M)
Export.FGT.Table(FGT_Results_Table = Acq_Ctx_2way.Early, TableName = "Umag.Ctx.Acq.Early", file.export = T)
# no significant effects

# Create Mid Trials mids object
A_Ctx_Mid <- subset(A_Ctx, trial %in% c(3,4)) 
unique(A_Ctx_Mid$trial)
A_Ctx_Mid_M <- means_EML_2way(A_Ctx_Mid)
# Analyze Mid Trials
Acq_Ctx_2way.Mid <- FGT_LMER_trialtype.Condition(dataset=A_Ctx_Mid_M)
Export.FGT.Table(FGT_Results_Table = Acq_Ctx_2way.Mid, TableName = "Umag.Ctx.Acq.Mid", file.export = T)
# no significant effects

# Create Late Trials mids object
A_Ctx_Late <- subset(A_Ctx, trial %in% c(5,6))
unique(A_Ctx_Late$trial)
A_Ctx_Late_M <- means_EML_2way(A_Ctx_Late)
# Analyze Late Trials
Acq_Ctx_2way.Late <- FGT_LMER_trialtype.Condition(dataset=A_Ctx_Late_M)
Export.FGT.Table(FGT_Results_Table = Acq_Ctx_2way.Late, TableName = "Umag.Ctx.Acq.Late", file.export = T)
# no significant effects

# No further follow-up analyses planned/needed


# ITI FPS (acquisition) --------------------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
A_ITI_Mean.log <-log.transform.mids(Umag.mids$Acq_ITI_Means)
Acq_ITI_LM <- FGT_LM.Condition(dataset=A_ITI_Mean.log)
Export.FGT.Table(FGT_Results_Table = Acq_ITI_LM, TableName = "Umag.ITI.Acq.Mean", file.export = T) 
# no significant effects

# Sensitivity analyses to check the effect of potential influential participants, identified via Assumption Checks: p 28, 107 (76 in most datasets)
A_ITI_Mean.log.noinf <- remove.influential.points(data_mids=A_ITI_Mean.log, influentialpoints=c(28,107,76))
Acq_ITI_noinf <- FGT_LM.Condition(dataset=A_ITI_Mean.log.noinf)
Export.FGT.Table(FGT_Results_Table = Acq_ITI_noinf, TableName = "Umag.ITI.Acq.noinf.Mean", file.export = T)
# Conclusion: significant change in result. Significant effect of condition following exclusion (28, 107, 76). Same for the exclusion of 28 and 107 alone. This discrepancy is described in the paper.

# Follow-up main effect condition:
Main_ITI_Acq<-with(expr=lm(FPS~1+Condition),data = A_ITI_Mean.log.noinf)
ACQ_ITI.followupCondition<- pool.means.condition(model=Main_ITI_Acq, m=100)
Export.FGT.FollowupCondition(FGT_FollowupResults=ACQ_ITI.followupCondition, 
                             TableName = "Umag.ITI.ACQ.FU.Condition", file.export = T)


# Day 2: Testing / Generalization phase -----------------------------------

# Habituation (Noise Alone) FPS (test) ---------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
T_Hab_Trials.log <- log.transform.mids(Umag.mids$Test_HAB_Trials)
#T_Hab_Trials.log<-Sensitivity_contingency(T_Hab_Trials.log) # Same results
Test_Hab_2way <- FGT_LMER_trial.Condition(dataset=T_Hab_Trials.log)
Export.FGT.Table(FGT_Results_Table = Test_Hab_2way, TableName = "Umag.Hab.Test.Trials", file.export = T)
# No significant effects (No follow-up analyses planned/needed)

# # Sensitivity analyses to check the effect of potential influential participants, identified via Assumption Checks: p24
T_Hab_Trials.log.noinf <- remove.influential.points(data_mids=T_Hab_Trials.log, influentialpoints=c(24))
Test_Hab_2way.noinf <- FGT_LMER_trial.Condition(dataset=T_Hab_Trials.log.noinf)
Export.FGT.Table(FGT_Results_Table = Test_Hab_2way.noinf, TableName = "Umag.Hab.Test.noinf", file.export = T)
# Conclusion: no significant change in result.


# Cued FPS Test -----------------------------------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
T_Cued_Trials.log <- log.transform.mids(Umag.mids$Test_Cued_Trials)
# T_Cued_Trials.log<-Sensitivity_contingency(T_Cued_Trials.log) # Same results
Test_Cued_3way <- FGT_LMER_Trialtype.Trial.Condition(dataset=T_Cued_Trials.log)
Export.FGT.Table(FGT_Results_Table = Test_Cued_3way, 
                 TableName = "Umag.Cue.Test.Trials", file.export = T) 
# Main effect Trialtype & Main effect Trial

# Sensitivity analyses to check the effect of potential influential participants, identified via Assumption Checks:  pp 26, 127 (and 100?)
T_Cued_Trials.log.noinf <- remove.influential.points(data_mids=T_Cued_Trials.log, influentialpoints=c(26, 127))
Test_Cued_3way.noinf <- FGT_LMER_Trialtype.Trial.Condition(dataset=T_Cued_Trials.log.noinf)
Export.FGT.Table(FGT_Results_Table = Test_Cued_3way.noinf, TableName = "Umag.Cue.Test.noinf", file.export = F)
# Conclusion: no significant change in result.

# Sensitivity analyses to check the influence of factor vs continuous variable "trailtype":
# Continuous analyses suggested by (http://www.frontiersin.org/Quantitative_Psychology_and_Measurement/10.3389/fpsyg.2015.00652/abstract) & https://www.sciencedirect.com/science/article/pii/S0005791618300612
T_Cued_Trials.log_Trialtypecontinuous<-Within.continuous.mids(T_Cued_Trials.log)
Test_Cued_3way_TTcontinuous <- FGT_LMER_Trialtype.Trial.Condition(dataset=T_Cued_Trials.log_Trialtypecontinuous)
Export.FGT.Table(FGT_Results_Table = Test_Cued_3way_TTcontinuous, TableName = "Umag.Cue.Test.Trials_TTcontinuous", file.export = F)
# Conclusion: no significant change in result.

# * Follow-up analyses 2: decompose factor "Trial" by analyzing means of early & late trials separately ####
T_Cue<-complete(T_Cued_Trials.log, action = "long", include = T) # Change mids to long format
unique(T_Cue$trial) # 6 trials

# Create Early Trials mids object
T_Cue_Early <- subset(T_Cue, trial %in% c(1,2)) 
unique(T_Cue_Early$trial)
T_Cue_Early_M <- means_EML_2way(T_Cue_Early)
# Analyze Early Trials
Test_Cue_2way.Early <- FGT_LMER_trialtype.Condition(dataset=T_Cue_Early_M)
Export.FGT.Table(FGT_Results_Table = Test_Cue_2way.Early, TableName = "Umag.Cue.Test.Early", file.export = T)
# Main effect Trialtype

# Follow-up Main effect trialtype (to determine direction of effect):
# Copy Relevant Model from LMER script:
T_CUE_Main_NO_2way.EARLY<-with(expr=lmer(FPS~1+Condition+trialtype+(1|pp), control=lmerControl(optimizer = 'bobyqa')),data = T_Cue_Early_M, REML=F)
# Call function to pool EMM:
Test_Cue_2way.Early.followupTrialtype <-pool.means.trialtype(model=T_CUE_Main_NO_2way.EARLY, m=100, phase="TEST") # Enter larger model from the specific comparison that is followed up. In this case Full main effects model.
Test_Cue_2way.Early.followupTrialtype # Note Qbar = pooled means (overall point estimate)
# Export results to table
Export.FGT.FollowupTrialtype(FGT_FollowupResults = Test_Cue_2way.Early.followupTrialtype, 
                             TableName = "Umag.Cue.Test.Early_FU.trialtype", file.export = T)
# The effects plotted
Trialtype.PlotD2(Test_Cue_2way.Early.followupTrialtype, rownames(Test_Cue_2way.Early.followupTrialtype), Test_Cue_2way.Early.followupTrialtype$Qbar)


# Create Mid Trials mids object
T_Cue_Mid <- subset(T_Cue, trial %in% c(3,4)) 
unique(T_Cue_Mid$trial)
T_Cue_Mid_M <- means_EML_2way(T_Cue_Mid)
# Analyze Mid Trials
Test_Cue_2way.Mid <- FGT_LMER_trialtype.Condition(dataset=T_Cue_Mid_M)
Export.FGT.Table(FGT_Results_Table = Test_Cue_2way.Mid, TableName = "Umag.Cue.Test.Mid", file.export = T)
# Main effect Trialtype

# Follow-up Main effect trialtype (to determine direction of effect):
# Copy Relevant Model from LMER script:
T_CUE_Main_NO_2way.MID<-with(expr=lmer(FPS~1+Condition+trialtype+(1|pp), control=lmerControl(optimizer = 'bobyqa')),data = T_Cue_Mid_M, REML=F)
# Call function to pool EMM:
Test_Cue_2way.MID.followupTrialtype <-pool.means.trialtype(model=T_CUE_Main_NO_2way.MID, m=100, phase="TEST") # Enter larger model from the specific comparison that is followed up. In this case Full main effects model.
Test_Cue_2way.MID.followupTrialtype # Note Qbar = pooled means (overall point estimate)
# Export results to table
Export.FGT.FollowupTrialtype(FGT_FollowupResults = Test_Cue_2way.MID.followupTrialtype, 
                             TableName = "Umag.Cue.Test.Mid_FU.trialtype", file.export = T)
# The effects plotted
Trialtype.PlotD2(Test_Cue_2way.MID.followupTrialtype, rownames(Test_Cue_2way.MID.followupTrialtype), Test_Cue_2way.MID.followupTrialtype$Qbar)


# Create Late Trials mids object
T_Cue_Late <- subset(T_Cue, trial %in% c(5,6)) 
unique(T_Cue_Late$trial)
T_Cue_Late_M <- means_EML_2way(T_Cue_Late)
# Analyze Late Trials
Test_Cue_2way.Late <- FGT_LMER_trialtype.Condition(dataset=T_Cue_Late_M)
Export.FGT.Table(FGT_Results_Table = Test_Cue_2way.Late, TableName = "Umag.Cue.Test.Late", file.export = T)
# Main effect Trialtype

# Follow-up Main effect trialtype (to determine direction of effect):
# Copy Relevant Model from LMER script:
T_CUE_Main_NO_2way.LATE<-with(expr=lmer(FPS~1+Condition+trialtype+(1|pp), control=lmerControl(optimizer = 'bobyqa')),data = T_Cue_Late_M, REML=F)
# Call function to pool EMM:
Test_Cue_2way.LATE.followupTrialtype <-pool.means.trialtype(model=T_CUE_Main_NO_2way.LATE, m=100, phase="TEST") # Enter larger model from the specific comparison that is followed up. In this case Full main effects model.
Test_Cue_2way.LATE.followupTrialtype # Note Qbar = pooled means (overall point estimate)
# Export results to table
Export.FGT.FollowupTrialtype(FGT_FollowupResults = Test_Cue_2way.LATE.followupTrialtype, 
                             TableName = "Umag.Cue.Test.Late_FU.trialtype", file.export = T)
# The effects plotted
Trialtype.PlotD2(Test_Cue_2way.LATE.followupTrialtype, rownames(Test_Cue_2way.LATE.followupTrialtype), Test_Cue_2way.LATE.followupTrialtype$Qbar)


# Contextual FPS Test -----------------------------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
T_Ctx_Trials.log <-log.transform.mids(Umag.mids$Test_Context_Trials)
#T_Ctx_Trials.log<-Sensitivity_contingency(T_Ctx_Trials.log) # Significant trial*condition interaction, no effect of Condition in follow-up analyses on epochs from early, mid and late trials (see below). No further follow-up.
Test_Ctx_3way <- FGT_LMER_Trialtype.Trial.Condition(dataset=T_Ctx_Trials.log)
Export.FGT.Table(FGT_Results_Table = Test_Ctx_3way, 
                 TableName = "Umag.Ctx.Test.Trials", file.export = T) 
# Main effect trial. No interactions (No follow-up analyses planned/needed)

# Sensitivity analyses to check the effect of potential influential participants, identified via Assumption Checks: pp 20, 31, 132
T_Ctx_Trials.log.noinf <- remove.influential.points(data_mids=T_Ctx_Trials.log, influentialpoints=c(20, 31, 132))
Test_Ctx_3way.noinf <- FGT_LMER_Trialtype.Trial.Condition.Workaround.D2.CTX(dataset=T_Ctx_Trials.log.noinf)
Export.FGT.Table(FGT_Results_Table = Test_Ctx_3way.noinf, TableName = "Umag.Ctx.Test.noinf", file.export = F)
# Conclusion: no significant change in result.

# Sensitivity analyses to check the influence of factor vs continuous variable "trailtype":
T_Ctx_Trials.log_Trialtypecontinuous<-Within.continuous.mids(T_Ctx_Trials.log)
Test_Ctx_3way_TTcontinuous <- FGT_LMER_Trialtype.Trial.Condition(dataset=T_Ctx_Trials.log_Trialtypecontinuous)
Export.FGT.Table(FGT_Results_Table = Test_Ctx_3way_TTcontinuous, TableName = "Umag.Ctx.Test.Trials_TTcontinuous", file.export = F)
# Conclusion: no significant change in result.

# * Follow-up analyses sensitivity analyses contingency: decompose factor "Trial" by analyzing means of early & late trials separately ####
T_Ctx<-complete(T_Ctx_Trials.log, action = "long", include = T) # Change mids to long format
unique(T_Ctx$trial) # 3 trials

# Create Early Trials mids object
T_Ctx_Early <- subset(T_Ctx, trial %in% c(1))
unique(T_Ctx_Early$trial)
T_Ctx_Early_M <- means_EML_2way(T_Ctx_Early)
# Analyze Early Trials
Test_Ctx_2way.Early <- FGT_LMER_trialtype.Condition(dataset=T_Ctx_Early_M)
Export.FGT.Table(FGT_Results_Table = Test_Ctx_2way.Early, TableName = "Umag.Ctx.Test.Early", file.export = T)
# No significant effects of condition for trial 1

# Create Mid Trials mids object
T_Ctx_Mid <- subset(T_Ctx, trial %in% c(2))
unique(T_Ctx_Mid$trial)
T_Ctx_Mid_M <- means_EML_2way(T_Ctx_Mid)
# Analyze Mid Trials
Test_Ctx_2way.Mid <- FGT_LMER_trialtype.Condition(dataset=T_Ctx_Mid_M)
Export.FGT.Table(FGT_Results_Table = Test_Ctx_2way.Mid, TableName = "Umag.Ctx.Test.Mid", file.export = T)
# No significant effects of condition for trial 2

# Create Late Trials mids object
T_Ctx_Late <- subset(T_Ctx, trial %in% c(3))
unique(T_Ctx_Late$trial)
T_Ctx_Late_M <- means_EML_2way(T_Ctx_Late)
# Analyze Late Trials
Test_Ctx_2way.Late <- FGT_LMER_trialtype.Condition(dataset=T_Ctx_Late_M)
Export.FGT.Table(FGT_Results_Table = Test_Ctx_2way.Late, TableName = "Umag.Ctx.Test.Late", file.export = T)
# No significant effects of condition for trial 3 (but main effect of trialtype here)

# Follow-up Main effect trialtype (to determine direction of effect):
# Copy Relevant Model from LMER script:
T_Ctx_Main_NO_2way.LATE<-with(expr=lmer(FPS~1+Condition+trialtype+(1|pp),control=lmerControl(optimizer="bobyqa")),data = T_Ctx_Late_M, REML=F)
# Call function to pool EMM:
Test_Ctx_2way.LATE.followupTrialtype <-pool.means.trialtype(model=T_Ctx_Main_NO_2way.LATE, m=100, phase="TEST") # Enter larger model from the specific comparison that is followed up. In this case Full main effects model.
Test_Ctx_2way.LATE.followupTrialtype # Note Qbar = pooled means (overall point estimate)
# Export results to table
Export.FGT.FollowupTrialtype(FGT_FollowupResults = Test_Ctx_2way.LATE.followupTrialtype,
                             TableName = "Umag.Ctx.Test.Late_FU.trialtype", file.export = T)
# The effects plotted
Trialtype.PlotD2(Test_Ctx_2way.LATE.followupTrialtype, rownames(Test_Ctx_2way.LATE.followupTrialtype), Test_Ctx_2way.LATE.followupTrialtype$Qbar)


# ITI FPS (Test) --------------------------------------------
# Note LMER assumptions checked after imputation (in 100 datasets) and satisfied after log(FPS)
T_ITI_Mean.log <-log.transform.mids(Umag.mids$Test_ITI_Means)
Test_ITI_LM <- FGT_LM.Condition(dataset=T_ITI_Mean.log)
Export.FGT.Table(FGT_Results_Table = Test_ITI_LM, TableName = "Umag.ITI.Test.Mean", file.export = F) 
# No significant effects

# Sensitivity analyses to check the effect of potential influential participants, identified via Assumption Checks: p 76, 107 (88 and 94 in most datasets)
T_ITI_Mean.log.noinf <- remove.influential.points(data_mids=T_ITI_Mean.log, influentialpoints=c(76, 107, 88, 94))
Test_ITI_LM_noinf <- FGT_LM.Condition(dataset=T_ITI_Mean.log.noinf)
Export.FGT.Table(FGT_Results_Table = Test_ITI_LM_noinf, TableName = "Umag.ITI.Test.noinf.Mean", file.export = T)
# Conclusion: no influence on results (no exclusion of pp)
