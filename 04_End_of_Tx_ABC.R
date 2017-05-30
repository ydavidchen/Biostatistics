################################################
# RACE data analysis: ABC Score Comparisons
# Primary endpoint = End of Treatment
# Script author: David (Youdinghuan) Chen
# Date: 05/19/2017
################################################

library(dplyr)
library(gdata)
library(geepack)
library(ggplot2)
library(grid)
library(gridExtra)
library(impute)
library(matrixStats)
library(reshape)

confint.geeglm <- function(object, parm, level=0.95, ...) {
  #'@description Method for `confint` on `geeglm` objects. Returns matrix with lower+upper CI.
  #'@param object geeglm using geepack package
  #'@param level Desired confidence level
  #'@author Ben Bolker via StackOverflow
  cc <- coef(summary(object));
  mult <- qnorm((1+level)/2);
  citab <- with(
    as.data.frame(cc),
    cbind(lwr=Estimate-mult*Std.err, upr=Estimate+mult*Std.err)
  );
  rownames(citab) <- rownames(cc);
  return(citab[parm, ])
}

setwd("~/Dropbox (Personal)/Dartmouth/ActiveStepRCT_data/")

#------------------------------------Dataset------------------------------------
## Load dataset:
RACE_data      <- read.csv("Race Data.csv")
RACE_data_dict <- read.xls("RACET Data Dictionary.xlsx")

## Endpoint: end of treatment:
RACE_data.end <- RACE_data[RACE_data$IntervalName=="End of Treatment", ]
RACE_data.end <- RACE_data.end[! is.na(RACE_data.end$abc_average), ]
nrow(RACE_data.end)

## Keep screened subjects with end-of-tx data only:
RACE_data.screen <- RACE_data[RACE_data$IntervalName=="Screening", ]
RACE_data.screen <- RACE_data.screen[RACE_data.screen$ReferenceID %in% RACE_data.end$ReferenceID, ]
RACE_data.screen <- RACE_data.screen[! is.na(RACE_data.screen$abc_average), ]
nrow(RACE_data.screen)

## Update post-screening data.frames:
RACE_data.end <- RACE_data.end[RACE_data.end$ReferenceID %in% RACE_data.screen$ReferenceID, ]
nrow(RACE_data.end)

identical(RACE_data.screen$ReferenceID, RACE_data.end$ReferenceID) #if FALSE, match first
identical(RACE_data.screen$ReferenceID, RACE_data.end$ReferenceID)

q1 <- ggplot(RACE_data.screen, aes(abc_average, fill=TreatmentCode, color=TreatmentCode)) +
  geom_density(alpha=0.25) + 
  scale_fill_manual(values=c("salmon", "darkolivegreen3")) +
  scale_color_manual(values=c("salmon", "darkolivegreen3")) +
  facet_wrap(~TreatmentCode) + 
  ggtitle("ABC score distribution at screening") + 
  theme_bw()

q3 <- ggplot(RACE_data.end, aes(abc_average, fill=TreatmentCode, color=TreatmentCode)) +
  geom_density(alpha=0.25) + 
  scale_fill_manual(values=c("salmon", "darkolivegreen3")) +
  scale_color_manual(values=c("salmon", "darkolivegreen3")) +
  facet_wrap(~TreatmentCode) + 
  ggtitle("ABC score distribution at End of Treatment") + 
  theme_bw()

grid.arrange(q1,q3,ncol=1)

#------------------------------------Pairwise t-tests------------------------------------
t.test(
  RACE_data.end$abc_average[RACE_data.end$TreatmentCode=="ACTIVESTEP"], 
  RACE_data.screen$abc_average[RACE_data.screen$TreatmentCode=="ACTIVESTEP"], 
  paired = TRUE
)

t.test(
  RACE_data.end$abc_average[RACE_data.end$TreatmentCode=="STANDARD"], 
  RACE_data.screen$abc_average[RACE_data.screen$TreatmentCode=="STANDARD"], 
  paired = TRUE
)

#------------------------------------GEE modeling------------------------------------
## Combine complete data from 2 time points:
combined_screen_end <- rbind(RACE_data.screen, RACE_data.end)
combined_screen_end <- combined_screen_end[order(combined_screen_end$ReferenceID),  ]
length(unique(combined_screen_end$ReferenceID)) #checkpoint: 350 expected

## Convert case-control status to 0,1:
combined_screen_end$TreatmentCode <- gsub("ACTIVESTEP", 1, combined_screen_end$TreatmentCode)
combined_screen_end$TreatmentCode <- gsub("STANDARD", 0, combined_screen_end$TreatmentCode)
combined_screen_end$TreatmentCode <- as.factor(combined_screen_end$TreatmentCode)

## Generate time-order variable:
combined_screen_end$timeorder <- ifelse(combined_screen_end$IntervalName=="Screening", 1, 2)
combined_screen_end$timeorder <- as.integer(combined_screen_end$timeorder)

gee2 <- geeglm(
  abc_average ~ TreatmentCode * timeorder, #note syntax
  data = combined_screen_end, 
  id = ReferenceID, 
  family = poisson("identity"), 
  corstr = "independence"
)
summary(gee2)$coef

gee2.tab <- as.data.frame(cbind(
  summary(gee2)$coef,
  confint(gee2)
))
colnames(gee2.tab)[colnames(gee2.tab) %in% c("lwr", "upr")] <- c("CI_lower", "CI.upper")
gee2.tab <- round(gee2.tab, 2)
gee2.tab$`Pr(>|W|)`[gee2.tab$`Pr(>|W|)`==0] <- "<2E-16"

grid.arrange(tableGrob(gee2.tab))
