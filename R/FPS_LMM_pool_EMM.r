# Pool LMER Estimated Marginal Means (EMM) of imputed data from the Fear Generalization Task (FGT) in the SAM study
# written by Rosalie Gorter, adapted by Milou Sep (function created)
# Last version 16.1.19

library(emmeans) # to calculate EMMs
library(ggplot2) # to visualize EMMs


# EMM Trialtype -----------------------------------------------------------

pool.means.trialtype <- function(model, m, phase){
  # model is the largest model, m is number of imputations, phase is "ACQ" or "TEST"
  
  results<-list()
  for(i in 1:m){#place emmeans() results for every imputation in list
    results[i]<-list(emmeans(model[[4]][[i]], 
                             pairwise ~ trialtype, #set contrasts
                             adjust="tukey",#p-value adjustment for multiple comparisons family wise (note: Benjamini-Hochberg correction is not appropriate here because p-values are not independent http://www-stat.wharton.upenn.edu/~steele/Courses/956/Resource/MultipleComparision/Writght92.pdf)
                             lmer.df = "kenward-roger"))
  }
  
  # # To check:
  # model[[4]][[1]]
  # results[[1]]
  # results[[1]]$emmeans
  # data.frame(results[[1]][1]$emmeans)$emmean 
  # data.frame(results[[i]][1]$emmeans)$SE
  
  # pool parameters according to rubins rules 
  # Rubin DB. Inference and missing data. Biometrika. 1976;63:581-92.
  # Marshall A, Altman DG, Holder RL, Royston P. Combining estimates of interest in prognostic modelling studies after multiple imputation: current practice and guidelines. BMC Med Res Methodol. 2009;9:57. doi:10.1186/1471-2288-9-57.
  # steps: calculate pooled parameter; calculate pooled standard error (pse); used pse to calculate 95% CI
  
  # Combine all imputed values into pooled mean and variance and 95% CI
  
  if (phase == "ACQ"){
  Qhat<-matrix(NA,ncol=m,nrow=2) #estimated means for every imputed dataset
  rownames(Qhat)<-c("mean.Threat","mean.Safe")
  colnames(Qhat)<-c(paste("imp.",1:m,sep=""))
  
  U_i<-matrix(NA,ncol=m,nrow=2) #estimated variance per imputed dataset associated with mean
  rownames(U_i)<-c("SE.Threat","SE.Safe")
  colnames(U_i)<-c(paste("imp.",1:m,sep=""))
  
  } else if (phase == "TEST"){
    Qhat<-matrix(NA,ncol=m,nrow=3) #estimated means for every imputed dataset
    rownames(Qhat)<-c("mean.Threat","mean.Safe", "mean.New")
    colnames(Qhat)<-c(paste("imp.",1:m,sep=""))
    
    U_i<-matrix(NA,ncol=m,nrow=3) #estimated variance per imputed dataset associated with mean
    rownames(U_i)<-c("SE.Threat","SE.Safe", "SE.New")
    colnames(U_i)<-c(paste("imp.",1:m,sep=""))
  }
  
  for(i in 1:m){ #place estimates per imputation in matrix
    Qhat[,i]<-data.frame(results[[i]][1]$emmeans)$emmean 
    U_i[,i]<-data.frame(results[[i]][1]$emmeans)$SE
  }
  
  #pool means and se's
  Ubar<-rowMeans(U_i) #within imputation variance
  B<-(1/(m-1)) * rowSums((Qhat-rowMeans(Qhat))^2)#between imputation variance
  Total.var <-Ubar+B #Total variance
  Qbar<-rowMeans(Qhat) #pooled means (overall point estimate)
  r_m<-(1+m^-1)*B/Ubar #relative increase in variance
  nu<-(m-1)*(1+r_m)^2 #degrees of freedom t reference distribution
  results[[1]]
  t.bound<-abs(qt(0.05/2, nu))#calculate t values for boundaries 95% CI based on adjusted df (nu)
  
  lower.95<-Qbar-t.bound*Total.var^(1/2) #rubin1987 p77
  upper.95<-Qbar+t.bound*Total.var^(1/2)
  
  pooled.results.means<-data.frame(cbind(Qbar,Total.var,lower.95,upper.95))
  return(pooled.results.means) #pooled estimates per Condition
}


# EMM Condition -----------------------------------------------------------

pool.means.condition <- function(model, m){
  # model is the largest model, m is number of imputations
  results<-list()
  for(i in 1:m){#place emmeans() results for every imputation in list
    results[i]<-list(emmeans(model[[4]][[i]], 
                             pairwise ~ Condition, #set contrasts
                             adjust="tukey",#p-value adjustment for multiple comparisons family wise (note: Benjamini-Hochberg correction is not appropriate here because p-values are not independent http://www-stat.wharton.upenn.edu/~steele/Courses/956/Resource/MultipleComparision/Writght92.pdf)
                             lmer.df = "kenward-roger"))
  }

  # pool parameters according to rubins rules 
  # Rubin DB. Inference and missing data. Biometrika. 1976;63:581-92.
  # Marshall A, Altman DG, Holder RL, Royston P. Combining estimates of interest in prognostic modelling studies after multiple imputation: current practice and guidelines. BMC Med Res Methodol. 2009;9:57. doi:10.1186/1471-2288-9-57.
  # steps: calculate pooled parameter; calculate pooled standard error; used pse to calculate 95% CI
  
  # Combine all imputed values into pooled mean and variance and 95% CI
  
  Qhat<-matrix(NA,ncol=m,nrow=3) #estimated means for every imputed dataset
  rownames(Qhat)<-c("mean.Delayed.Stress","mean.Direct.Stress", "mean.No.Stress")
  colnames(Qhat)<-c(paste("imp.",1:m,sep=""))
  
  U_i<-matrix(NA,ncol=m,nrow=3) #estimated variance per imputed dataset associated with mean
  rownames(U_i)<-c("SE.Delayed.Stress","SE.Direct.Stress", "SE.No.Stress")
  colnames(U_i)<-c(paste("imp.",1:m,sep=""))
  
  for(i in 1:m){ #place estimates per imputation in matrix
    Qhat[,i]<-data.frame(results[[i]][1]$emmeans)$emmean 
    U_i[,i]<-data.frame(results[[i]][1]$emmeans)$SE
  }
  
  #pool means and se's
  Ubar<-rowMeans(U_i) #within imputation variance
  B<-(1/(m-1)) * rowSums((Qhat-rowMeans(Qhat))^2)#between imputation variance
  Total.var <-Ubar+B #Total variance
  Qbar<-rowMeans(Qhat) #pooled means (overall point estimate)
  r_m<-(1+m^-1)*B/Ubar #relative increase in variance
  nu<-(m-1)*(1+r_m)^2 #degrees of freedom t reference distribution
  results[[1]]
  t.bound<-abs(qt(0.05/2, nu))#calculate t values for bounderies 95% CI based on adjusted df (nu)
  
  lower.95<-Qbar-t.bound*Total.var^(1/2) #rubin1987 p77
  upper.95<-Qbar+t.bound*Total.var^(1/2)
  
  pooled.results.means<-data.frame(cbind(Qbar,Total.var,lower.95,upper.95))
  return(pooled.results.means) #pooled estimates per Condition
}


# Function to plot mean epochs follow-up analyses -------------------------
Trialtype.PlotD2 <- function(dataset_to_plot, x, y){
  plot2<-
    ggplot(data=dataset_to_plot, aes(x, y)) +
    theme_classic() +
    ylab('FPS \n (mean+/-95%CI)') + xlab("") +   # axis labels & title
    scale_x_discrete(limits=c("mean.Threat","mean.Safe","mean.New"), # Change order x-axis
                     labels=c("Threat ", "Safe", "New ")) +
    scale_fill_manual( # Colours: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
      values= c( Safe = "navyblue", Threat = "red3", New = "black")) +
    # labels =c(Delayed = "delayed-stress",
    #           Direct = "direct-stress",
    #           Control = "no-stress")) +
    # adjust fontsize legend
    theme(legend.title=element_text(size=12), legend.text=element_text(size=9))+
    # Add mean & CI's
    geom_bar(stat="identity", position = position_dodge(.9), alpha=.8) +
    geom_errorbar(aes(ymin=lower.95, ymax=upper.95), width=.2, position=position_dodge(.9))
  
  return(plot2)
}