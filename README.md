# SAM_FGT
Data Imputation and Analyses of Fear Generalization Task in SAM study, as described in:
- Sep, M.S.C., Gorter, R., van Ast, V.A., Joëls, M., & Geuze, E. (2019) No Time-Dependent Effects of Psychosocial Stress on Fear Contextualization and Generalization: A Randomized-Controlled Study With Healthy Participants. Chronic Stress, 3, 247054701989654. https://doi.org/10.1177/2470547019896547

#### Step 1: Data
Datasets (available on Dataverse):
- `SAM_FGT.csv` contains fear-potentiated startle responses (FPS) during the FGT task
- `SAM_FGT_Amperage.csv` contains shock intensities used in the conditioning phase of the FGT, per participant.
- `SAM_questionnaires.sav` contains questionnaire information (including fear contingency scores)
- `SAM_Codes_Task_Protocol_Versions.csv` contains information on the task versions, per participant

#### Step 2: Descriptive Statistics
`FGT_descritpvies.R` (in the 'R' folder) loads:
- `SAM_FGT.csv` and `SAM_FGT_Amperage.csv` for the description of shock intensities and missing values.
- `SAM_questionnaires.sav` and `SAM_Codes_Task_Protocol_Versions.csv` to prepare data for sensitivity analyses with fear contingency scores (saved as `participants.for.contingency.sensitivity.analyses.rda` in the folder `processed_data`)

#### Step 3: Multiple Imputation (MI) of FPS
The `FPS_imputation.R` script (in the 'R' folder) loads `SAM_FGT.csv` and prepares the data for multiple imputation (the cleaned data is saved as `INPUT_IMPUTATIE_FGT_Umag.rda` in the folder `processed_data`).

Perform MI via the function `Impute_FGT_EMG_SAM` to:
  1) to **impute individual trials** (set `sorttype` to 'trials'). 
     - Note: trials will be imputed *if more than 1/3 of the trials in a category is present* (in other words if <2/3 missing), if less than 1/3 is present (in other words if >2/3 is missing; `missing code 4`) all the trials (for that category) will be set to missing 
  2) to **create imputed means** (set `sorttype` to 'mean').
     - Note: means will be based on imputed trials if *more than 2/3 of trials is present* (in other words if <1/3 missing; `missing code 1`), or imputed directly if less than 2/3 of the trials is present (or in other words, *if >1/3 missing*; `missing code 2` ).

The **imputed datasets** are saved in the folder `processed_data` as:
  1) `OUPUT_IMPUTATIE_FGT_out_M50_MAXIT100_Umag_AllMeans.rda`
  2) `OUPUT_IMPUTATIE_FGT_out_M50_MAXIT100_Umag_Trials.rda`

Note, the function `Impute_FGT_EMG_SAM` also saves the output of the script automatically with a generic name -`OUPUT_IMPUTATIE_FGT_out_endofscript.rda`- in the folder `processed_data`.

#### Step 4: MI quality checks & data transformations
- `FPS_mids_Quality.R` (in the 'R' folder) creates a sorted `mids` object -`Umag.mids.28.01.19.rda`- from the imputed data in the folder `processed_data`, that is used for the analyses.

#### Step 5: Assumptions linear mixed effect models (LMM):
The assumptions for LMM analyses are checked within each imputed dataset via `FPS_LMM_Assumptions_call.R`  (in the 'R' folder). Note, this script loads `Umag.mids.28.01.19.rda` and renders the script `FPS_LMM_Assumptions_source.R` to a rmarkdown file for each outcome measure (files will appear in 'R' folder).

#### Step 6: LMM analyses:
The LMM analyses are performed within each imputed dataset via `FPS_LMM_Analyses_call.R`  (in the 'R' folder). This script loads `Umag.mids.28.01.19.rda` and sources:

1) scripts with **analyses functions**:
    - `FPS_LMM_Trialtype.Trial.Condition.r`
    - `FPS_LMM_trialtype.Condition.r`
    - `FPS_LMM_trial.Condition.R`
    - `FPS_LMM_LM_Condition.r`

2) a script to pool (and plot) **LMM estimates**: `FPS_LMM_pool_EMM.r`
3) a script to **transform `mids` objects**: `FGT_mids_transformations.R`
4) a script to **export results**: `FPS_LMM_Results_to_Word.R`

#### Step 7: Data Visualization:
- Figures: `FPS_LMM_Results_to_Plot.R`
