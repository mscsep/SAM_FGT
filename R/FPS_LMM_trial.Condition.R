# Supporting function for the analysis of the data from the Fear Generalization Task (FGT) in the SAM study
# Written by Milou Sep

FGT_LMER_trial.Condition <- function(dataset){

  {
    ## 2way interactions ##
    # Define Models
    FullModel <-with(expr=lmer(FPS~1+Condition*trial+(1|pp),control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F) 
    Main_NO_2way<-with(expr=lmer(FPS~1+Condition+trial+(1|pp),control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F)
    # Compare Models
    Condition.trial <- pool.compare(FullModel, Main_NO_2way, method='wald')
    
    ## Main effects ##
    # Define Models
    Main_Condition<-with(expr=lmer(FPS~1+trial+(1|pp),control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F)
    Main_trial<-with(expr=lmer(FPS~1+Condition+(1|pp),control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F)
    # Compare Models
    Condition <- pool.compare(Main_NO_2way,Main_Condition, method='wald')
    trial <- pool.compare(Main_NO_2way,Main_trial, method='wald')
  }
  
  # Table of model comparisons (Publication) --------------------------------
  twoway.C.T <- c("Condition x Trial", Condition.trial$Dm, Condition.trial$rm, Condition.trial$df1, Condition.trial$df2, Condition.trial$pvalue)
  Main.effect.Condition <- c("Condition", Condition$Dm, Condition$rm, Condition$df1, Condition$df2, Condition$pvalue)
  Main.effect.trial <- c("Trial", trial$Dm, trial$rm, trial$df1, trial$df2,  trial$pvalue)
  
  WORD.table.FGT <- data.frame(rbind(twoway.C.T ,
                                     Main.effect.Condition, Main.effect.trial))
  colnames(WORD.table.FGT) <- c("Model", "Dm", "rm", "df1", "df2", "pvalue")
  
  # Data cells need to be changed to numeric (not factor), otherwise problems with flextable().  
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.character) # First convert to character ..
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.numeric) # .. then to numeric
  
  # Return function output --------------------------------------------------
  
  return(WORD.table.FGT)
  
}
