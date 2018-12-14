# Phase II Clinical Trial Analysis: PFS & OS Curves in 1 Figure
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Date: 06/05/18
# Notes:

rm(list=ls());
library(ggsci);
library(survival);
library(survminer);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions_folbrite.R");

#--------------------------------------------Survival Objects--------------------------------------------
DAYS_PER_MONTH <- 30.44; 

## Load data:
data <- load_data();

## Median time to last follow-up for the entire study:
followup.median.os <- median(data$Overall.Survival.Days..n.30., na.rm=TRUE) / DAYS_PER_MONTH; 
followup.median.os
sum(is.na(data$Overall.Survival.Days..n.30.))

followup.median.pfs <- median(data$PFS.Days..n.36., na.rm=TRUE) / DAYS_PER_MONTH; 
followup.median.pfs
sum(is.na(data$PFS.Days..n.36.))

## Curve objects:
pfsDat <- subset(data, ! is.na(PFS.Days..n.36.) );
pfsDat$years <- pfsDat$PFS.Days..n.36. / 365.25;
pfsDat$pfs.status <- as.logical(pfsDat$Event.Free.vs..Progression..n.39.);
pfs <- survfit(Surv(years, pfs.status) ~ 1, data=pfsDat); 
pfs

osDat <- subset(data, ! (is.na(Overall.Survival.Days..n.30.) | is.na(Overall.Survival.Events..n.30.)) );
osDat$years <- osDat$Overall.Survival.Days..n.30. / 365.25;
osDat$os.status <- osDat$Overall.Survival.Events..n.30. == 3; 
os <- survfit( Surv(years, os.status) ~ 1, data=osDat );
os

#--------------------------------------------Statement in manuscript--------------------------------------------
## Answer to Erick's question 1a:
print(os, print.rmean=TRUE)
summary(os, times=followup.median.os / 12)
# summary(os, times=48 / 12)

## Answer to Erick's question 1b:
print(pfs, print.rmean=TRUE)
summary(pfs, times=followup.median.pfs / 12)
# summary(pfs, times=48 / 12)

#--------------------------------------------Data visualization--------------------------------------------
fitsCombined <- list(PFS=pfs, OS=os); 

mySurvTheme$legend.title <- element_blank();
ggsurvplot(
  fitsCombined,
  combine = TRUE,
  risk.table = FALSE,
  tables.theme = theme_cleantable(),  #clean risk table
  conf.int = FALSE,
  conf.int.style = "step", #"step" or "ribbon"
  censor = TRUE, # emove censor points
  size = 1.5,
  alpha = 0.75, #rel. position
  xlab = "Time (years)",
  legend = c(0.8, 0.2),
  palette = "jco",
  linetype = "strata",
  ggtheme = mySurvTheme
) 

ggsave("~/Downloads/Overall_OS_and_PFS.png", dpi=300, width=8, height=8)

