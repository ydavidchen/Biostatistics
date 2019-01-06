# Phase II Clinical Trial Analysis: PFS & OS Curves in One Figure
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Copyright (c) 2018 @ydavidchen
# Dates: 06/05/18; 01/06/19 (revision)
# Notes:

rm(list=ls());
if(interactive()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions_folbrite.R");

## Load data:
data <- load_data();

#--------------------------------------------Survival Objects--------------------------------------------
## Median time to last follow-up for the entire study:
lastFollowUp.median <- median(data$followUpDurationInDays) / DAYS_PER_MONTH;
lastFollowUp.median 

## Curve objects:
pfs <- survfit(Surv(pfs_years, pfs_status) ~ 1, data=data); 
pfs

os <- survfit(Surv(os_years, os_status) ~ 1, data=data);
os

#--------------------------------------------Statement in manuscript--------------------------------------------
## Answer to Erick's question 1a:
print(os, print.rmean=TRUE)
summary(os, times=lastFollowUp.median/12)
summary(os, times=48/12)

## Answer to Erick's question 1b:
print(pfs, print.rmean=TRUE)
summary(pfs, times=lastFollowUp.median/12)
summary(pfs, times=48/12)

#--------------------------------------------Data visualization--------------------------------------------
fitsCombined <- list(PFS=pfs, OS=os);
MY_SURV_THEME$legend.title <- element_blank();

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
  ggtheme = MY_SURV_THEME
) 

ggsave("~/Downloads/Lansigan_et_al_Figures/Combined_OS_and_PFS.png", dpi=300, width=8, height=8)
