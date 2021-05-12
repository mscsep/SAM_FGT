# Supporting function for the analysis of the data from the Fear Generalization Task (FGT) in the SAM study
# Written by Milou Sep

FGT_LM.Condition <- function(dataset){
  
  # Compare Conditions
  Base<-with(expr=lm(FPS~1),data = dataset)
  Main_Condition<-with(expr=lm(FPS~1+Condition),data = dataset)
  # Compare models
  Condition <- pool.compare(Main_Condition, Base, method='wald')
  
  # Make vectors for data in each row
  Main.effect.Condition <- c("Condition", Condition$Dm, Condition$rm, Condition$df1, Condition$df2, Condition$pvalue)
  
  WORD.table.FGT <- data.frame(rbind( Main.effect.Condition))
  colnames(WORD.table.FGT) <- c("Model", "Dm", "rm", "df1", "df2", "pvalue")
  
  # Data cells need to be changed to numeric (not factor), otherwise problems with flextable().  
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.character) # First convert to character ..
  WORD.table.FGT[c(2,3,4,5,6)] <- sapply(WORD.table.FGT[c(2,3,4,5,6)], as.numeric) # .. then to numeric

  # Return function output --------------------------------------------------
  
  return(WORD.table.FGT)
  
}
