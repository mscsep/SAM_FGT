# Supporting function to export the results from pool.compare() in the analysis of the data from the Fear Generalization Task (FGT) in the SAM study to word tables 
# written by Milou Sep
# version 13.9.18

# Load required packages --------------------------------------------------
library(flextable) # To make publication ready tables for word | info: https://davidgohel.github.io/flextable/articles/overview.html
library(officer) # To export chi tables to word | info: https://davidgohel.github.io/flextable/articles/layout.html   

# Function for Chi2 Tables ------------------------------------------------
Export.FGT.Table <- function(FGT_Results_Table, TableName, file.export){
  # Create Flextable
  FGT_Table <- flextable(FGT_Results_Table)
  # Add Layout settings to flextable().
  # Make significant value's Bold [Info: https://cran.r-project.org/web/packages/flextable/vignettes/format.html#bold]
  FGT_Table <-  bold(FGT_Table, i= ~pvalue < 0.05 , j="pvalue")  # Note is j= is included, only value not total row bold
  # set number of digits  # https://davidgohel.github.io/flextable/articles/format.html#set_formatter-function
  FGT_Table <-  display(FGT_Table, col_key="df1", pattern = "{{df1}}", formatters = list(df1~sprintf("%.00f",df1))) # no digits for df (note if 1 digit would be required, use %.01f etc.).
  FGT_Table <-   theme_vanilla(FGT_Table)  # remove thick lines
  # FGT_Table <-  set_header_labels(FGT_Table, LogLikelihood="Log Likelihood" ,deltadf= "delta df", pvalue="p-value") # change names
  FGT_Table <-   autofit(FGT_Table) # adapt table dimensions to content
  
  if (file.export == T){
    doc <- read_docx()
    doc <- body_add_flextable(doc, value = FGT_Table, align="center")
    print(doc, target = paste0("results/",TableName, date(),".docx"))
  }
  
  return(FGT_Table)
}


Export.FGT.FollowupTrialtype <- function(FGT_FollowupResults, TableName, file.export){
  names<-rownames(FGT_FollowupResults)
  # Create Flextable
  FGT_Table <- flextable(data=cbind(names, FGT_FollowupResults))
  FGT_Table<-set_header_labels(FGT_Table, names="Trialtype", Qbar = "Pooled Means", Total.var = "variance",
                               lower.95="lower 95CI", upper.95 = "upper 95CI") # Note Qbar = pooled means (overall point estimate)
  # Add Layout settings to flextable().
  FGT_Table <-   theme_vanilla(FGT_Table)
  FGT_Table <-   autofit(FGT_Table)
  
  if (file.export == T){
    doc <- read_docx()
    doc <- body_add_flextable(doc, value = FGT_Table, align="center")
    print(doc, target = paste0("results/", TableName, date(),".docx"))
  }
  
  return(FGT_Table)
}


Export.FGT.FollowupCondition <- function(FGT_FollowupResults, TableName, file.export){
  names<-rownames(FGT_FollowupResults)
  # Create Flextable
  FGT_Table <- flextable(data=cbind(names, FGT_FollowupResults))
  FGT_Table<-set_header_labels(FGT_Table, names="Condition", Qbar = "Pooled Means", Total.var = "variance",
                               lower.95="lower 95CI", upper.95 = "upper 95CI")
  # Add Layout settings to flextable ().
  FGT_Table <-   theme_vanilla(FGT_Table)
  FGT_Table <-   autofit(FGT_Table)
  
  if (file.export == T){
    doc <- read_docx()
    doc <- body_add_flextable(doc, value = FGT_Table, align="center")
    print(doc, target = paste0("results/", TableName, date(),".docx"))
  }
  
  return(FGT_Table)
}
