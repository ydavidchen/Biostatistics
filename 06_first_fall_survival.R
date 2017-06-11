#############################################################################################
# RACE RCT data  analysis: Time intervals till the first fall
# David Chen
# Start date: 05/28/2017
# Notes:
#############################################################################################

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

#------------------------------------Extract binary fall/no fall for each time point------------------------------------
## Survival analysis with censoring:
RACE_data.surv <- RACE_data[ , c("ReferenceID", "TreatmentCode","IntervalName","fallen_past_3months")]
RACE_data.surv <- RACE_data.surv[! RACE_data.surv$IntervalName %in% c("End of Treatment","Repeat Treatment"), ]
RACE_data.surv$fallen_past_3months[is.na(RACE_data.surv$fallen_past_3months)] <- 0

## Code in time point in numeric format for sorting:
RACE_data.surv$timepoint <- NA;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="Screening"         ] <-  0;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="3-Month Follow-up" ] <-  3;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="6-Month Follow-up" ] <-  6;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="9-Month Follow-up" ] <-  9;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="Annual Follow-up"  ] <- 12;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="15-Month Follow-up"] <- 15;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="18-Month Follow-up"] <- 18;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="21-Month Follow-up"] <- 21;
RACE_data.surv$timepoint[RACE_data.surv$IntervalName=="2-Year Follow-up"  ] <- 24

RACE_data.surv <- RACE_data.surv[with(RACE_data.surv, order(ReferenceID, timepoint)), ] # Sort by vector name [ReferenceID] then [timepoint]
RACE_data.surv <- droplevels(RACE_data.surv)

activeStepIDs <- unique(RACE_data.surv$ReferenceID[RACE_data.surv$TreatmentCode=="ACTIVESTEP"]) #for query later

## For rows with 1 for a given patient, drop rows after the first event:
last_timePt <- "(no event)"; #for people never fell
patients <- unique(RACE_data.surv$ReferenceID);
first_fall_rec <- NULL;
for(id in patients){
  ## Select a subject's records:
  subj_indices <- which(RACE_data.surv$ReferenceID %in% id)
  DFi <-  RACE_data.surv[subj_indices, ]
  
  ## Conditional update: if a patient ever fell...
  if(c(1) %in% DFi$fallen_past_3months){
    IntervalName <- DFi$IntervalName[DFi$fallen_past_3months==1][1]
  }else{
    IntervalName <- last_timePt
  }
  
  ## Update results matrix:
  first_fall_rec <- rbind(
    first_fall_rec, 
    cbind(id, IntervalName)
  )
}

first_fall_rec <- as.data.frame(first_fall_rec);

first_fall_rec$timepoint <- NA;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="Screening"         ] <-  0;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="3-Month Follow-up" ] <-  3;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="6-Month Follow-up" ] <-  6;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="9-Month Follow-up" ] <-  9;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="Annual Follow-up"  ] <- 12;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="15-Month Follow-up"] <- 15;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="18-Month Follow-up"] <- 18;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="21-Month Follow-up"] <- 21;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="2-Year Follow-up"  ] <- 24;
first_fall_rec$timepoint[first_fall_rec$IntervalName=="(no event)"        ] <- 24 #or Inf?

## Code in events:
first_fall_rec$event <- ifelse(
  test = first_fall_rec$IntervalName == "(no event)", 
  yes  = FALSE,
  no   = TRUE
)
first_fall_rec$TreatmentCode <- ifelse(
  test = first_fall_rec$id %in% activeStepIDs,
  yes  = "ACTIVESTEP",
  no   = "STANDARD"
)
colnames(first_fall_rec)

## Kaplan-Meier Curve:
with(first_fall_rec, Surv(timepoint, event))
mySurv <- Surv(first_fall_rec$timepoint, event=first_fall_rec$event)
myFit  <- survfit(mySurv ~ TreatmentCode, data=first_fall_rec)
summary(myFit)
survdiff(mySurv ~ TreatmentCode, data=first_fall_rec) #log-rank test
# autoplot(myFit, conf.int=FALSE) + 
#   labs(x="Time in months", y="Survival probability", title="Comparison of time-to-event") +
#   scale_color_brewer(palette="Set1") +
#   theme_bw()
ggsurvplot(
  myFit, 
  data = first_fall_rec, 
  risk.table = TRUE,
  palette = "Set1", 
  pval = TRUE, 
  xlab = "Time in months"
)

