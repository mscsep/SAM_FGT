#' ---	
#' title: "Visualization of Raw FPS response during the Fear Generalization Task (FGT) in the SAM study"
#' author: "Milou Sep"	
#' date: august 2019
#' output:	
#'   html_document: default	
#' ---	

# clear workspace
rm(list=ls())

# packages ----------------------------------------------------------------
# library(readr) # to load .txt data
library(dplyr)
library(tidyr)
library(reshape2) # to change wide to long format
library(Rmisc) # to calculate mean and CI in raw data
library(ggplot2) # to create plots
library(ggpubr) # To arrange plots.

# load and restructure data -----------------------------------------------
FGT_batch <- read.csv2("data/SAM_FGT.csv", na.strings = c("NaN","5555","8888","9999"))
FGT_batch$subjName <- gsub("SMGH","", FGT_batch$subjName)
Umag <- subset.data.frame(FGT_batch, select = c(grep("Umag", names(FGT_batch)),Condition, subjName))

# Change name "Condition" in "Group"
Umag%>%dplyr::rename(Group=Condition) -> Umag

# 1) Change data to long format
melt(Umag, id.vars = c("subjName", "Group"), measure.vars = grep("Umag", names(Umag)),
     variable.name = 'variable', value.name='FPS') ->Umag_long
str(Umag_long)

# 2) Split variable names
separate(Umag_long, col =variable, into = c("phase", "probe", NA, "trial"), sep = "_", convert = T) -> Umag_long
str(Umag_long) 
unique(Umag_long$phase)

# 3) Split trialtype code (letter) and trial number (digit)
Umag_long %>% mutate(trialtype = gsub("[[:digit:]]","",trial), # Note = replace letter by ""
                     trial = as.numeric(gsub("[^[:digit:]]","",trial))) -> Umag_long  # Note = replace not a letter by ""

# 4) Add HAB and ITI as trialtypes
Umag_long %>% mutate(
  trialtype = ifelse(probe == "ITI", "ITI", trialtype),
  trialtype = ifelse(probe == "HAB", "HAB", trialtype),
  trialtype = ifelse(probe == "HABctx", "HABctx", trialtype),
) -> Umag_long #  %>% filter(probe == "Cue")  # to check
# to check
Umag_long[which(Umag_long$probe == 'ITI'),] # Check iti
Umag_long[ which( Umag_long$trialtype == "HABctx"),]  # for HAB!

# 5) Change to factor
cols <- c("subjName","Group","phase", "probe", "trialtype")
#Use lapply() to coerce and replace the chosen columns:
Umag_long[cols] <- lapply(Umag_long[cols], factor)
# checking
unique(Umag_long$probe) #-> Cue, Ctx, ITI, Hab, Habctx
unique(Umag_long$trialtype) 

# Add labels to Group and trialtype for readability in plots (and change order of appearance in plots)
Umag_long$probe <- factor(Umag_long$probe, levels=c("HAB", "HABctx", "Cue", "ITI", "Ctx"),
                          labels = c("NAh", "pre-acq ctx", "CUE", "ITI", "CTX"))
attributes(Umag_long$probe)

attributes(Umag_long$Group)
Umag_long$Group <- factor(Umag_long$Group, levels=c(3,2,1),
                          labels =c("No-Stress", "Immediate-Stress", "Delayed-Stress"))

attributes(Umag_long$trialtype)
Umag_long$trialtype <- factor(Umag_long $trialtype, levels=c( "T" ,"S" ,  "N"  , "ITI", "HAB", "HABctx"  ),
                              labels = c("CTX+", "CTX-", "G-CTX", "ITI", "NAh", "pre-acq context"))
# 6) Change to numeric
Umag_long$FPS<-as.numeric(Umag_long$FPS)

# Calculate mean and CI ---------------------------------------------------
group.CI(FPS~Group+phase+probe+trial+trialtype, Umag_long, ci = 0.95) -> means_long
str(means_long)

# Function to plot raw data -----------------------------------------------
Plot_Raw_FGT <- function(data){
  plot<- 
    ggplot(data = data, aes(x=trial, y= FPS.mean, group = factor(trialtype), colour= trialtype)) +
    # Layout: 
    facet_grid(.~ Group ) + # Groups next to each other
    theme_classic() +
    scale_x_continuous(breaks = unique (data$trial)) +
    xlab("Trialnumber") +
    ylab('FPS (mean +/-95%CI)') +  # ax labels & title
    labs(colour = "Trialtype")+
    # adjust font size legend and axis
    theme(legend.title=element_text(size=12), legend.text=element_text(size=9),
          axis.text = element_text(size=9), axis.title= element_text(size=10))+
    #  add color to plot: https://stackoverflow.com/questions/17180115/manually-setting-group-colors-for-ggplot2
    scale_color_manual( # Colours: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
      values= c('NAh' = 'black',
                'CTX+' = 'red2',
                'CTX-' = 'mediumblue',
                'G-CTX' = 'forestgreen',
                'ITI' = 'gray50' )) +
    # add Data: 
    geom_line(aes(colour=trialtype), stat="summary" , fun.y="mean", na.rm=T, size=1) + # remove.na=TRUE is used to remove warnings about missing values (which are present in raw data)
    geom_point(stat="summary",fun.y="mean", na.rm=T, size=2)  + # Add datapoints
    # .. and errorbars (info: https://stackoverflow.com/questions/48800212/set-error-bars-to-standard-deviation-on-a-ggplot2-bar-graph)
    geom_errorbar(aes(ymin=FPS.lower, ymax=FPS.upper), position=position_dodge(.5), width=0.5)   
  
  return(plot)}


# Call function to plot data per phase ------------------------------------

# Acquisition: Noise Alone (Figure 4A)
means_long %>%
  filter(phase =='sA',
         probe == 'NAh') %>% 
  Plot_Raw_FGT() %>%
  + ggtitle(("Pre-aquisition noise alone trials")) ->acq_NAh

# CUE (Figure 4B)
means_long %>%
  mutate(trial = ifelse(probe == "ITI", (trial*2), trial)) %>%  # Change the trial numbers of the ITI trials (i.e. *2) for visualization, so that it is clear that ITI probes were delivered throughout acquisition phase (not just the first part)
  filter(phase =='sA',
         probe == 'CUE'|
           probe == 'ITI') %>% 
  Plot_Raw_FGT() %>%
  + ggtitle(("Contextualization of cued fear")) -> acq_cue

# CTX (Figure 4C)
means_long %>%
  filter(phase =='sA',
         probe == 'CTX'|
           probe == 'ITI') %>% 
  Plot_Raw_FGT() %>%
  + ggtitle(("Contextual fear expression (acquisition)"))  -> acq_ctx


# arrange the plots of the acquisition phase in Figure 4
FGT_acq.plots <- 
  ggarrange(
    acq_NAh,
    acq_cue,
    acq_ctx,
    labels = c("A", "B", "C"),
    ncol = 1, nrow = 3,
    font.label=list(size=12, face="bold"),
    align = "v",
    legend="right",
    common.legend = F)

annotate_figure(FGT_acq.plots, 
                fig.lab=c("Figure 4"), 
                fig.lab.pos = "top.right",
                fig.lab.face = "bold",
                fig.lab.size = 12)

ggsave(paste0("results/Figure4.Acq.Plots.colour.",date(),".eps"), device="eps", dpi = 800, height = 7, width = 7, limitsize = T )



# TEST: Noise Alone, Cue, CTX ---------------------------------------------

# NAh; Figure 5A
means_long %>%
  filter(phase =='sG',
         probe == 'NAh') %>% 
  Plot_Raw_FGT() %>%
  + ggtitle(("Pre-memory test noise alone trials")) ->test_NAh

# CUE; Figure 5B
means_long %>%
  mutate(trial = ifelse(probe == "ITI", (trial*2), trial)) %>%  # Change the trial numbers of the ITI trials (i.e. *2) for visualization, so that it is clear that ITI probes were delivered throughout acquisition phase (not just the first part)
  filter(phase =='sG',
         probe == 'CUE'|
           probe == 'ITI') %>% 
  Plot_Raw_FGT() %>%
  + ggtitle(("Context-dependency of cued fear memory")) -> test_cue

# CTX; Figure 5C
means_long %>%
  filter(phase =='sG',
         probe == 'CTX'|
           probe == 'ITI') %>% 
  Plot_Raw_FGT() %>%
  + ggtitle(("Contextual fear expression (memory)"))  -> test_ctx

# arrange the plots of the test phase in Figure 5.
FGT_test.plots <- 
  ggarrange(
    test_NAh,
    test_cue,
    test_ctx,
    labels = c("A", "B", "C"),
    ncol = 1, nrow = 3,
    font.label=list(size=12, face="bold"),
    align = "v",
    legend="right",
    common.legend = F)

annotate_figure(FGT_test.plots, 
                fig.lab=c("Figure 5"), 
                fig.lab.pos = "top.right",
                fig.lab.face = "bold",
                fig.lab.size = 12)

ggsave(paste0("results/Figure5.Test.Plots.colour.",date(),".eps"), device="eps", dpi = 800, height = 7, width = 7, limitsize = T )
