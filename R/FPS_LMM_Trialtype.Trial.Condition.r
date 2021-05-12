# Supporting function for the analysis of the data from the Fear Generalization Task (FGT) in the SAM study
# Written by Milou Sep

FGT_LMER_Trialtype.Trial.Condition <- function(dataset){
  
  # Note, the Bobyaqa optimizer method was used for better model convergence (Note, this optimizer was also used for MCT analyses (https://github.com/mscsep/SAM_MCT).
  # Info on LMER convergence problems: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
  # Info on optimizer checks http://svmiller.com/blog/2018/06/mixed-effects-models-optimizer-checks/
  # ?allFit() # Note: allFit() function doesn't work with Mids objects.
  
  FullModel <-with(expr=lmer(FPS~1+Condition*trialtype*trial+(1|pp),control=lmerControl(optimizer="bobyqa") ),data = dataset, REML=F) 
  
  Full_NO_3way<-with(expr=lmer(FPS~1+Condition+trialtype+trial+ # Main effects
                                 Condition:trialtype+Condition:trial+trialtype:trial+ # 2 way interaction effects
                                 (1|pp) ,control=lmerControl(optimizer="bobyqa")),data = dataset,REML=F)  # Random effects term
  
  ### 2 Way interactions [compared to Full Model without 3Way] ###
  # Define Models
  Full_NO_Condition.trial<-with(expr=lmer(FPS~1+Condition+trialtype+trial+ 
                                            Condition:trialtype+trialtype:trial+
                                            (1|pp) ,control=lmerControl(optimizer="bobyqa") ),data = dataset, REML=F)
  
  Full_NO_trial.trialtype<-with(expr=lmer(FPS~1+Condition+trialtype+trial+
                                            Condition:trialtype+Condition:trial+
                                            (1|pp) ,control=lmerControl(optimizer="bobyqa") ),data = dataset, REML=F)
  
  Full_NO_Condition.trialtype<-with(expr=lmer(FPS~1+Condition+trialtype+trial+
                                                Condition:trial+trialtype:trial+
                                                (1|pp),control=lmerControl(optimizer="bobyqa") ),data = dataset, REML=F)
  
  # 3 way interaction
  Condition.trial.trialtype <-  pool.compare(FullModel,Full_NO_3way, method='wald') # To test if there is a 3way interaction, and the Log-Likelihood (and df) of the full model  # pool.compare(main, base, method='wald')
  
  # 2 way interactions
  Condition.trial <- pool.compare(Full_NO_3way,Full_NO_Condition.trial, method='wald')
  trial.trialtype <- pool.compare(Full_NO_3way,Full_NO_trial.trialtype, method='wald')
  Condition.trialtype <-  pool.compare(FullModel, Full_NO_Condition.trialtype, method='wald')
  
  ### Main effects ###
  # Define Models  
  Main_NO_2ways<-with(expr=lmer(FPS~1+Condition+trialtype+trial+(1|pp)  ,control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F)
  Main_Condition<-with(expr=lmer(FPS~1+trial+trialtype+(1|pp) ,control=lmerControl(optimizer="bobyqa") ),data = dataset, REML=F)
  Main_trial<-with(expr=lmer(FPS~1+Condition+trialtype+(1|pp) ,control=lmerControl(optimizer="bobyqa") ),data = dataset, REML=F)
  Main_trialtype<-with(expr=lmer(FPS~1+Condition+trial+(1|pp) ,control=lmerControl(optimizer="bobyqa") ),data = dataset, REML=F)
  
  # Compare Models
  Condition <- pool.compare(Main_NO_2ways,Main_Condition, method='wald')
  trial <- pool.compare(Main_NO_2ways,Main_trial, method='wald')
  trialtype <- pool.compare(Main_NO_2ways,Main_trialtype, method='wald')
  
  
  # Table of model comparisons (for publication) --------------------------------
  # Outcomes of pool.compare() (RG: export " Dm		rm		df1	df2		p " to table)
  #?pool.compare
  
  # Make vectors for data in each row
  threeway.C.T.T <- c("Condition x Trial x Trialtype", Condition.trial.trialtype$Dm, Condition.trial.trialtype$rm, Condition.trial.trialtype$df1, Condition.trial.trialtype$df2, Condition.trial.trialtype$pvalue)
  
  twoway.C.T <- c("Condition x Trial", Condition.trial$Dm, Condition.trial$rm, Condition.trial$df1, Condition.trial$df2, Condition.trial$pvalue)
  twoway.C.Tt <- c("Condition x Trialtype", Condition.trialtype$Dm, Condition.trialtype$rm, Condition.trialtype$df1, Condition.trialtype$df2, Condition.trialtype$pvalue)
  twoway.T.Tt <- c("Trial x Trialtype", trial.trialtype$Dm, trial.trialtype$rm, trial.trialtype$df1, trial.trialtype$df2, trial.trialtype$pvalue)
  
  Main.effect.Condition <- c("Condition", Condition$Dm, Condition$rm, Condition$df1, Condition$df2, Condition$pvalue)
  Main.effect.trialtype <- c("Trialtype", trialtype$Dm, trialtype$rm, trialtype$df1, trialtype$df2,  trialtype$pvalue)
  Main.effect.trial <- c("Trial", trial$Dm, trial$rm, trial$df1, trial$df2,  trial$pvalue)
  
  WORD.table.FGT <- data.frame(rbind(threeway.C.T.T ,
                                     twoway.C.T, twoway.C.Tt , twoway.T.Tt,
                                     Main.effect.Condition, Main.effect.trialtype, Main.effect.trial))
  colnames(WORD.table.FGT) <- c("Model", "Dm", "rm", "df1", "df2", "pvalue")
  
  # Data cells need to be changed to numeric (not factor), otherwise problems with flextable().
  # str(WORD.table.FGT)
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.character) # First convert to character ..
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.numeric) # .. then to numeric
  # sapply(WORD.table.FGT, class) # check data class of table.
  
  # Return function output --------------------------------------------------
  
  return(WORD.table.FGT)
  
}