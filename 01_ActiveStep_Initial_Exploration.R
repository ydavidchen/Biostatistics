####################################################
# Data exploration: ActiveStep RCT
# Script author: David (Youdinghuan) Chen
# Date: 04/25/2017
####################################################

rm(list=ls())
library(dplyr)
library(gdata)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(reshape)

# Load data set:
RACE_data      <- read.csv("~/Dropbox (Personal)/Dartmouth/ActiveStepRCT_data/Race Data.csv")
RACE_data_dict <- read.xls("~/Dropbox (Personal)/Dartmouth/ActiveStepRCT_data/RACET Data Dictionary.xlsx")

# View variables
table(RACE_data_dict$Variable.Group.Name)
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Follow-up:  Fall Details"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Demographics"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Follow-up:  Falls"]

# Number of unique individuals:
patients <- unique(RACE_data$ReferenceID)

# Count number of visit of each unique participant:
# num_visits <- c()
# for(p in patients){
#   s <- sum(RACE_data$ReferenceID == as.character(p)) 
#   num_visits  <- c(num_visits, s)
# }
# names(num_visits) <- patients
# hist(num_visits)

# table(RACE_data$SiteName)
# table(RACE_data$IntervalName)

## Identify and check (unique) Subjects underwent screening
# RACE_screen <- RACE_data$ReferenceID[RACE_data$IntervalName=="Screening"] 
# anyDuplicated(RACE_screen) #checkpoint: 0 required
# all(RACE_screen %in% patients) #checkpoint: TRUE required
# RACE_endTx <- RACE_data$ReferenceID[RACE_data$IntervalName=="End of Treatment"]
# all(RACE_endTx %in% RACE_screen) #checkpoint: TRUE required

## Set blank entries to NA:
RACE_data$discontinuation_date[RACE_data$discontinuation_date==""] <- NA


nrow(RACE_data)
RACE_data <- RACE_data[is.na(RACE_data$discontinuation_date), ]
nrow(RACE_data)
table(RACE_data$IntervalName)
# Cases are subjects who are assigned to the ActiveStep treatment group:
# The primary endpoint is 3 month period, which does NOT have abc_average
RACE_cases <- RACE_data[(RACE_data$TreatmentCode=="ACTIVESTEP") & (RACE_data$IntervalName=="End of Treatment"), ]
RACE_ctrls <- RACE_data[(RACE_data$TreatmentCode=="STANDARD")   & (RACE_data$IntervalName=="End of Treatment"), ]

anyDuplicated(RACE_cases$ReferenceID) #checkpoint tp ensure no duplicates: 0 expected
anyDuplicated(RACE_ctrls$ReferenceID) #checkpoint to ensure no duplicates: 0 expected

# Alternative hypothesis 1: At the end of the study period, ActiveStep group has higher ABC average score than controls
wilcox.test(
  RACE_cases$abc_average,
  RACE_ctrls$abc_average,
  alternative = "greater",
  correct = TRUE
)

# Alternative hypothesis: At the 3-month follow-up, ActiveStep group has lower number of Falls
RACE_cases <- RACE_data[(RACE_data$TreatmentCode=="ACTIVESTEP") & (RACE_data$IntervalName=="3-Month Follow-up"), ]
RACE_ctrls <- RACE_data[(RACE_data$TreatmentCode=="STANDARD")   & (RACE_data$IntervalName=="3-Month Follow-up"), ]
anyDuplicated(RACE_cases$ReferenceID) #checkpoint tp ensure no duplicates: 0 expected
wilcox.test(
  RACE_cases$times_fallen_3_months,
  RACE_ctrls$times_fallen_3_months,
  alternative = "two.sided",
  correct = TRUE
)

RACE_data1 <- RACE_data[RACE_data$IntervalName=="3-Month Follow-up", ]
anyDuplicated(RACE_data1 $ReferenceID) #check duplicates

# Histograms of times fallen:
times_fall <- RACE_data1  %>%
  select(ReferenceID, TreatmentCode, times_fallen_3_months) %>%
  group_by(TreatmentCode, tier=cut(times_fallen_3_months, breaks=c(1,3,Inf), include.lowest=TRUE)) %>% 
  summarise(n=n())
times_fall <- times_fall[!is.na(times_fall$tier), ]
ggplot(times_fall , aes(x=tier, y=n, fill=TreatmentCode)) +
  geom_bar(stat="identity", position="dodge") +  
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Total times fallen", y="Number of participants", title="Total number of falls at 3-month follow-up") +
  theme_bw()

# Boxplot by case-control status
sum(!is.na(RACE_data1$times_fallen_3_months[RACE_data1$TreatmentCode=="ACTIVESTEP"]))
sum(!is.na(RACE_data1$times_fallen_3_months[RACE_data1$TreatmentCode=="STANDARD"]))

ggplot(RACE_data1 , aes(TreatmentCode, times_fallen_3_months, color=TreatmentCode)) +
  geom_boxplot() + 
  scale_color_brewer(palette="Set1") +
  labs(y="Times fallen", title="Primary endpoint: Number of times fallen at 3-month follow-up \n (n=138 of 422 patients with data)") +
  theme_bw()

RACE_data2 <- RACE_data[RACE_data$IntervalName=="Screening", ]
# Check distribution of age: 
qplot(c(RACE_cases$age_at_screening, RACE_ctrls$age_at_screening),geom="histogram", bins=15) + 
  geom_vline(xintercept=75, col="red", linetype="dashed") + 
  labs(title="Age at screening (all unique patients, both arms combined)", x="age") + 
  theme_bw()

ggplot(RACE_data2, aes(TreatmentCode, age_at_screening, color=TreatmentCode)) + 
  geom_boxplot() + 
  scale_color_brewer(palette="Set1") +
  labs(y="Age at screening", title="Distribution of ABC Scores at screening") +
  theme_bw()

# Histogram stratified by age: using 75y/o as a cutoff: 
age_at_screen <- RACE_data2 %>%
  select(ReferenceID, TreatmentCode, age_at_screening) %>%
  group_by(TreatmentCode, tier=cut(age_at_screening, breaks=c(0,75,Inf), include.lowest=TRUE)) %>% 
  summarise(n=n())

ggplot(age_at_screen, aes(x=tier, y=n, fill=TreatmentCode)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Age group (cut-off = 75 y)", y="Number of participants", title="Age at screening") +
  theme_bw()

peripheral_neural <- RACE_data2 %>%
  select(ReferenceID, TreatmentCode, peripheral_neuropathy) %>%
  group_by(TreatmentCode, peripheral_neuropathy) %>% 
  summarise(n=n())
peripheral_neural$peripheral_neuropathy <- gsub(1, "yes", peripheral_neural$peripheral_neuropathy)
peripheral_neural$peripheral_neuropathy <- gsub(0, "no", peripheral_neural$peripheral_neuropathy)

# Histogram of number falls caused injury by group:
ggplot(peripheral_neural, aes(x=peripheral_neuropathy, y=n, fill=TreatmentCode)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Presence (1=yes, 0=no)", y="Number of participants", title="Presence of peripheral neuropathy, \n at the beginning (i.e.screening)") +
  theme_bw()

# Number of falls leading to injury by group, at the time of screening
ggplot(RACE_data2, aes(TreatmentCode,peripheral_neuropathy)) + 
  geom_bar(stat="identity") + 
  labs(y="Number of individuals", title="Number of Individuals with peripheral neuropathy at screening") +
  scale_y_continuous(expand=c(0, 0)) + #strip white space
  theme_bw()

# Select time point of interest:
RACE_data2 <- RACE_data[RACE_data$IntervalName=="Screening", ]

# Histogram by tiers:
stroke_presence <- RACE_data2 %>%
  select(ReferenceID, TreatmentCode, stroke) %>%
  group_by(TreatmentCode, stroke) %>% 
  summarise(n=n())
stroke_presence$stroke <- gsub(1, "yes", stroke_presence$stroke)
stroke_presence$stroke <- gsub(0, "no", stroke_presence$stroke)
ggplot(stroke_presence, aes(x=stroke, y=n, fill=TreatmentCode)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Presence", y="Number of participants", title="Presence of Stroke, \n at the beginning (i.e.screening)") +
  theme_bw()

# Number of individuals with stroke by case-control status:
ggplot(RACE_data2, aes(TreatmentCode, stroke)) + 
  geom_bar(stat="identity") + 
  labs(y="Number of individuals", title="Number of Individuals with Stroke at screening") +
  scale_y_continuous(expand=c(0, 0)) + #strip white space
  theme_bw()

# Select time point of interest:
RACE_data2 <- RACE_data[RACE_data$IntervalName=="Screening", ]

# Histogram by tiers: 
abc_baseline <- RACE_data2 %>%
  select(ReferenceID, TreatmentCode, abc_average) %>%
  group_by(TreatmentCode, tier=cut(abc_average, breaks=c(1,63, Inf), include.lowest=TRUE)) %>% 
  summarise(n=n())
abc_baseline <- abc_baseline[! is.na(abc_baseline$tier), ]
ggplot(abc_baseline, aes(x=tier, y=n, fill=TreatmentCode)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(expand=c(0, 0)) +
  labs(x="Group (median ABC average = 63)", y="Number of participants", title="Average Activity-specific Balance Confidence (ABC) scale, \n at the beginning (i.e.screening)") +
  theme_bw()

# Distribution of ABC scores
qplot(RACE_data2$abc_average, geom="histogram") +
  labs(x="ABC Average", title="Average Activity-specific Balance Confidence (ABC) at screening") + 
  theme_bw()

# ABC score by case-control status, at baseline:
ggplot(RACE_data2, aes(abc_average, fill=TreatmentCode)) + 
  geom_histogram(bins=15) + 
  facet_wrap(~TreatmentCode) +
  labs(x="ABC Score", title="Distribution of ABC Scores: Histograms at baseline") +
  scale_x_discrete(expand=c(0, 0)) + #strip white space
  scale_y_continuous(expand=c(0, 0))

ggplot(RACE_data2, aes(TreatmentCode, abc_average, color=TreatmentCode)) + 
  geom_boxplot() + 
  scale_color_brewer(palette="Set1") +
  labs(y="Mean ABC Score", title="Distribution of ABC Scores at baseline") +
  theme_bw()

# Covariates to keep in mind:
RACE_data3 <- RACE_data[RACE_data$IntervalName=="End of Treatment", ]
ggplot(RACE_data3, aes(TreatmentCode, total_num_sessions, color=TreatmentCode)) + 
  geom_boxplot() + 
  scale_color_brewer(palette="Set1") +
  labs(y="Number of sessions", title="Total number of sessions used at the end of study period") +
  theme_bw()