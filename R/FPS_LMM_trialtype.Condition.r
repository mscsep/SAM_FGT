# Supporting function for the analysis of the data from the Fear Generalization Task (FGT) in the SAM study
# Written by Milou Sep

FGT_LMER_trialtype.Condition <- function(dataset){

  {
    ## 2way interactions ##
    # Define Models
    FullModel <-with(expr=lmer(FPS~1+Condition*trialtype+(1|pp), control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F) 
    Main_NO_2way<-with(expr=lmer(FPS~1+Condition+trialtype+(1|pp),control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F)
    # Compare Models
    Condition.trialtype <- pool.compare(FullModel, Main_NO_2way, method='wald')
    
    ## Main effects ##
    # Define Models
    Main_Condition<-with(expr=lmer(FPS~1+trialtype+(1|pp),control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F)
    Main_trialtype<-with(expr=lmer(FPS~1+Condition+(1|pp),control=lmerControl(optimizer="bobyqa")),data = dataset, REML=F)
    # Compare Models
    Condition <- pool.compare(Main_NO_2way,Main_Condition, method='wald')
    trialtype <- pool.compare(Main_NO_2way,Main_trialtype, method='wald')
  }
  
  # Table of model comparisons (Publication) --------------------------------
  twoway.C.T <- c("Condition x Trialtype", Condition.trialtype$Dm, Condition.trialtype$rm, Condition.trialtype$df1, Condition.trialtype$df2, Condition.trialtype$pvalue)
  Main.effect.Condition <- c("Condition", Condition$Dm, Condition$rm, Condition$df1, Condition$df2, Condition$pvalue)
  Main.effect.trialtype <- c("Trialtype", trialtype$Dm, trialtype$rm, trialtype$df1, trialtype$df2,  trialtype$pvalue)
  
  WORD.table.FGT <- data.frame(rbind(twoway.C.T ,
                                     Main.effect.Condition, Main.effect.trialtype))
  colnames(WORD.table.FGT) <- c("Model", "Dm", "rm", "df1", "df2", "pvalue")
  
  # Data cells need to be changed to numeric (not factor), otherwise problems with flextable().  
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.character) # First convert to character ..
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.numeric) # .. then to numeric
  
  # Return function output --------------------------------------------------
  
  return(WORD.table.FGT)
  
}
