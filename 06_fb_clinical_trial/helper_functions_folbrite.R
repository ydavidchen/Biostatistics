# Phase II Clinical Trial Analysis: helper functions
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Copyright (c) 2018 @ydavidchen
# Date: (ongoing)
# Notes:

library(ggplot2);
library(ggsci);
library(survival);
library(survminer);
library(doParallel); registerDoParallel(detectCores() - 1);

DAYS_PER_YEAR <- 365.25;
DAYS_PER_MONTH <- 30.44;

load_data <- function(path="~/Dropbox (Christensen Lab)/Christensen Lab - 2018/StatsConsulting_2018/Survival_curve_proj/Sheet1_Updated.csv") {
  #'@description Load latest survival data and pre-process for data analysis
  #'@param path Path to latest data, resaved as CSV file
  ## Basic data loading
  df <- read.csv(path, stringsAsFactors=FALSE);
  df[df == ""] <- NA; 
  stopifnot(! any(duplicated(df$Patient.ID)));
  
  df$DATE.OF.ENROLLMENT <- as.Date(df$DATE.OF.ENROLLMENT, "%m/%d/%y"); 
  df$Last.Date.of.FU <- as.Date(df$Last.Date.of.FU, "%m/%d/%y");
  df$followUpDurationInDays <- as.numeric(df$Last.Date.of.FU - df$DATE.OF.ENROLLMENT);
  
  ## Separate data preprocessing for OS and PFS
  df$os_status[df$Overall.Survival.Events..n.30. == 3] <- 1; #here, death event is 3, alive is 2
  df$os_status[df$Overall.Survival.Events..n.30. == 2] <- 0;
  df$os_status[is.na(df$os_status)] <- 0; #assign censored status

  df$pfs_status <- df$Event.Free.vs..Progression..n.39.;
  
  ## Assign time to censored observations:
  df$os_days_complete <- df$Overall.Survival.Days..n.30.;
  df$os_days_complete[is.na(df$os_days_complete)] <- df$followUpDurationInDays[is.na(df$os_days_complete)]; 
  df$os_years <- df$os_days_complete / DAYS_PER_YEAR;
  
  df$pfs_days_complete <- df$PFS.Days..n.36.;
  df$pfs_days_complete[is.na(df$pfs_days_complete)] <- df$followUpDurationInDays[is.na(df$pfs_days_complete)]; 
  df$pfs_years <- df$pfs_days_complete / DAYS_PER_YEAR;
  
  return(df);
}

sentenceCase <- function(x) {
  #'@description Helper function to convert a string to sentence case, useful for graphic titles
  #'@param x A character string
  substr(x, 1, 1) <- toupper(substr(x, 1, 1));
  return(x);
}

## Graphic themes:
MY_SURV_THEME <- theme_classic2() +
    theme(axis.text=element_text(size=20,color="black"), 
          axis.title=element_text(size=20,color="black"),
          legend.title=element_text(size=16,color="black",face="bold"), legend.text=element_text(size=16,color="black",face="bold"), 
          legend.box.background=element_rect(colour="black",size=2),
          text=element_text(size=20),
          strip.text.x=element_text(size=15,colour="black",face="bold"));
