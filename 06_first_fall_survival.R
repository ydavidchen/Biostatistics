#############################################################################################
# RACE RCT data  analysis: Time intervals till the first fall
# David Chen
# Start date: 05/28/2017
# Notes:
# 1. Disregard everything after the first fall or first no-show
# 2. However, for subjects with both event (first fall) & no-show: prioritize the event
#############################################################################################

rm(list=ls())

library(dplyr)
library(gdata)
library(ggplot2)
library(ggfortify)
library(grid)
library(gridExtra)
library(matrixStats)
library(reshape)
library(survival)
library(survminer) #CRAN

setwd("~/Dropbox (Personal)/Dartmouth/ActiveStepRCT_data/")

#--------------------------------------------Dataset--------------------------------------------
## Load dataset:
RACE_data      <- read.csv("Race Data.csv", stringsAsFactors=FALSE)
RACE_data_dict <- read.xls("RACET Data Dictionary.xlsx")

## Check variables in dictionary:
table(RACE_data_dict$Variable.Group.Name)
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Follow-up:  Fall Details"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Demographics"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Follow-up:  Falls"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Baseline Assessment"]


RACE_data.screen <- RACE_data[RACE_data$IntervalName=="Screening", ]
RACE_data.rest <- RACE_data[RACE_data$IntervalName!="Screening", ]
mean(unique(RACE_data.rest$ReferenceID) %in% RACE_data.screen$ReferenceID)

#------------------------------------Extract binary fall/no fall for each time point------------------------------------
## Survival analysis with censoring:
RACE_data.surv <- RACE_data[ , c("ReferenceID", "TreatmentCode","IntervalName","fallen_past_3months")]
RACE_data.surv <- RACE_data.surv[! RACE_data.surv$IntervalName %in% c("End of Treatment","Repeat Treatment"), ]
RACE_data.surv$fallen_past_3months[is.na(RACE_data.surv$fallen_past_3months)] <- 0
patients <- unique(RACE_data.surv$ReferenceID); #all unique patients
activeStepIDs <- unique(RACE_data.surv$ReferenceID[RACE_data.surv$TreatmentCode=="ACTIVESTEP"]) #for query later

## Code in time point in numeric format for sorting:
RACE_data.surv$IntervalName.num <- NA;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="Screening"         ] <-  0;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="3-Month Follow-up" ] <-  3;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="6-Month Follow-up" ] <-  6;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="9-Month Follow-up" ] <-  9;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="Annual Follow-up"  ] <- 12;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="15-Month Follow-up"] <- 15;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="18-Month Follow-up"] <- 18;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="21-Month Follow-up"] <- 21;
RACE_data.surv$IntervalName.num[RACE_data.surv$IntervalName=="2-Year Follow-up"  ] <- 24;
RACE_data.surv$IntervalName.num <- as.numeric(RACE_data.surv$IntervalName.num);

RACE_data.surv <- RACE_data.surv[with(RACE_data.surv, order(ReferenceID, IntervalName.num)), ]; # Sort by vector name [ReferenceID] then [IntervalName.num]
RACE_data.surv <- droplevels(RACE_data.surv);

RACE_data.screen <- RACE_data.surv[RACE_data.surv$IntervalName=="Screening", ]
RACE_data.rest <- RACE_data.surv[RACE_data.surv$IntervalName != "Screening", ]
mean(RACE_data.rest$ReferenceID %in% RACE_data.screen$ReferenceID)

## Add in censoring time:
RACE_data.surv$CensoringTime <- NA; 
for(id in patients){
  df <- RACE_data.surv[RACE_data.surv$ReferenceID==id, ];
  RACE_data.surv$CensoringTime[RACE_data.surv$ReferenceID==id] <- df$IntervalName.num[nrow(df)] #last available record
}
class(RACE_data.surv$CensoringTime) #"numeric" required

## Data frame for survival analysis:
last_timePt <- Inf; #identifier (of type numeric) for people never fell
first_fall_rec <- NULL;
for(id in patients){
  ## Select a subject's records:
  subj_indices <- which(RACE_data.surv$ReferenceID %in% id);
  DFi <-  RACE_data.surv[subj_indices, ];
  
  ## Conditional update: if a patient ever fell...
  ## The first instance of fall will be the "event encounter"
  if(c(1) %in% DFi$fallen_past_3months){
    EventTime <- as.numeric(DFi$IntervalName.num[DFi$fallen_past_3months == 1][1]);
  }else{
    EventTime <- last_timePt;
  }
  
  ## Conditional update: if a patient ever failed to show up:
  ## The first instance of missing will be censored time
  CensoringTime <- as.numeric(unique(DFi$CensoringTime)); # should only have 1 censoring time
  
  ## Update results matrix:
  first_fall_rec <- rbind.data.frame(
    first_fall_rec, 
    cbind.data.frame(id, EventTime, CensoringTime)
  );
};

## Add in group:
first_fall_rec$TreatmentCode <- ifelse(
  test = first_fall_rec$id %in% activeStepIDs,
  yes  = "ACTIVESTEP",
  no   = "STANDARD"
); 

## Define events by comparing Censoring vs. Event time
first_fall_rec$event <- NA;
for(i in 1:nrow(first_fall_rec)){
  ## First code in events solely based on EventTime
  if(first_fall_rec$EventTime[i] == Inf){
    first_fall_rec$event[i] <- 0; #"alive"
  } else {
    first_fall_rec$event[i] <- 1; #"dead"
  };
  
  ## Second, if a patient had both 1+ fall or no-show:
  ## 1) Whether censored or not
  ##    24mo is the "administrative censoring" and does not count!
  ##     0mo is the exception!
  ## 2) If both an Event & censoring occurred during the study period, Event should take priority over censoring
  ## We are interested in those cases that are NOT administratively censored
  if( #If during this study period an event ever occurred AND censoring ever occurred
    ( ! first_fall_rec$CensoringTime[i] %in%  c(0,24) ) & 
    ( first_fall_rec$EventTime != Inf                 )
  ){
    first_fall_rec$event[i] <- 1; #"dead"
  } else if ( #If censoring occurred, but no fall, then considered alive
    ( ! first_fall_rec$CensoringTime[i] %in%  c(0,24) ) & 
    ( first_fall_rec$EventTime == Inf                 )
  ){
    first_fall_rec$event[i] <- 0; #"dead"
  };
};
colnames(first_fall_rec);
# first_fall_rec$start <- 0; # needed for interval censoring
class(first_fall_rec$start);
class((first_fall_rec$EventTime));
class((first_fall_rec$event))

## Kaplan-Meier Curve:
## Use default settings rather than interval
rm(EventTime) #avoid confusion
first_fall_rec$EventTime[first_fall_rec$EventTime==Inf] <- 24; #restore to enable plotting
mySurv <- Surv(first_fall_rec$EventTime, event=first_fall_rec$event);
mySurv
myFit  <- survfit(mySurv ~ TreatmentCode, data=first_fall_rec);
summary(myFit)
survdiff(mySurv ~ TreatmentCode, data=first_fall_rec) #log-rank test

autoplot(myFit, conf.int=FALSE, censor=FALSE) +
  labs(x="Time in months", y="Survival probability", title="Comparison of time-to-event") +
  scale_color_brewer(palette="Set1") +
  theme_bw()

ggsurvplot(
  myFit, 
  data = first_fall_rec, 
  risk.table = TRUE, #show risk table
  censor = FALSE, #don't draw bars
  pval = TRUE,  #show log-rank test P-value
  xlab = "Time in months",
  palette = "Set1"
)
