# Descriptive statistics of the data from the Fear Generalization Task (FGT) in the SAM study
# Written by Milou Sep.

# load packages
library(dplyr)
library(Rmisc) # For Mean and CI calculation of raw data
library(haven) # to load SPSS file in r

# Shock Intensities -------------------------------------------------------
# Load FGT data
FGT_batch <- read.csv2("data/SAM_FGT.csv", na.strings = c("NaN","5555","8888","9999"))
# load shock amperage values
FGT_SWU <- read.csv2("data/SAM_FGT_Amperage.csv")

# Select only participants that completed the FGT (n=117)
FGT_participants <- FGT_batch$subjName
# Changes characters to patterns for grepl(). Info on matching multiple patterns https://stackoverflow.com/questions/6947587/matching-multiple-patterns
FGT_pattern<-paste(FGT_participants, collapse ="|") 
FGT_SWU <- FGT_SWU[grepl(FGT_pattern, as.character(FGT_SWU$Participantnummer)),]

# Calculate mean, standard deviation, range of shock intensities
summary(FGT_SWU$`Stroomsterkte..mA.`) # mean & range
sd(FGT_SWU$`Stroomsterkte..mA.`) # standard deviation
# Calculate means CI / Condition
SWU.Mean.Condition.CI <- group.CI(`Stroomsterkte..mA.`~Conditie, FGT_SWU, ci = 0.95)
oneway.test(`Stroomsterkte..mA.`~Conditie, FGT_SWU) # NS


# Count missing FPS data per type missing (for paper) ----------------------------
# Note Missing codes (provided with Matlab code):
# 5555 % Digital registration error, message on screen & MissingValue 5555.
# 8888 % Excessive baseline activity. (More than 2 SD deviation)
# 9999 % Latency onset not valid (Note, this only affects onset scores, so not present in magnitude variables)

FGT_batch_count.missing <- read.csv2("data/SAM_FGT.csv")
Umag_for_missing <- subset.data.frame(FGT_batch_count.missing, select = c(grep("Umag", names(FGT_batch_count.missing))))
# Number of missing values in original data
n.missing.Technical<-sum(Umag_for_missing == '5555') #= 798
n.missing.Noise<-sum(Umag_for_missing == '8888') #= 256
# Total Number of observations
n.observations<-nrow(Umag_for_missing)*ncol(Umag_for_missing)
# percent.missing
(n.missing.Technical/n.observations)*100  # = 7.105
(n.missing.Noise/n.observations)*100  # = 2.279 %
# 0-responses
n.nulltesponses<-sum(Umag_for_missing == "0") #= 240
(n.nulltesponses/n.observations)*100 # 2.137 %


# FGT Contingency ---------------------------------------------------------
# load data
SAM_questionnaires <- read_sav("data/SAM_questionnaires.sav") # Questionnaires
SAM_versions <- read.csv("data/SAM_Codes_Task_Protocol_Versions.csv") # Information on task versions
# select required variables
SAM_questionnaires %>% select(SubjectNumber,Condition,shock_indicator)->SAM_questionnaires2 # NOTE, 59 missing values in 'shock indicator', because this question was added to the questionnaire later. 
SAM_versions %>% select(SubjectNumber,FGTversion)->SAM_versions2
# merge information
full_join(SAM_versions2, SAM_questionnaires2, by="SubjectNumber")->FGT_contingency
# Create variable to indicate if threat context was identified correctly.
FGT_contingency %>% mutate(Correct_FGT_contingency = case_when(
  shock_indicator == FGTversion ~ 1, # correct response
  shock_indicator != FGTversion ~ 0 # incorrect response
))-> FGT_contingency
# Differences between experimental groups?
group.CI(Correct_FGT_contingency~Condition, FGT_contingency, ci = 0.95)
oneway.test(Correct_FGT_contingency~Condition, FGT_contingency) # NS
# save participants with correct responses
FGT_contingency%>% filter(Correct_FGT_contingency == 1) %>% select(SubjectNumber) %>% saveRDS(.,"processed_data/participants.for.contingency.sensitivity.analyses.rda")
