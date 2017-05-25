################################################
# RACE data analysis: ABC Score Comparisons
# Primary endpoint = Annual Follow-up
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

## Annual endpoint:
RACE_data.ann <- RACE_data[RACE_data$IntervalName=="Annual Follow-up", ]
RACE_data.ann <- RACE_data.ann[! is.na(RACE_data.ann$abc_average), ]
nrow(RACE_data.ann)

## Keep screened subjects with end-of-tx data only:
RACE_data.screen <- RACE_data[RACE_data$IntervalName=="Screening", ]
RACE_data.screen <- RACE_data.screen[RACE_data.screen$ReferenceID %in% RACE_data.ann$ReferenceID, ]
RACE_data.screen <- RACE_data.screen[! is.na(RACE_data.screen$abc_average), ]
nrow(RACE_data.screen)

## Update post-screening data.frames:
RACE_data.ann <- RACE_data.ann[RACE_data.ann$ReferenceID %in% RACE_data.screen$ReferenceID, ]
nrow(RACE_data.ann)

identical(RACE_data.screen$ReferenceID, RACE_data.ann$ReferenceID) #if FALSE, match first
identical(RACE_data.screen$ReferenceID, RACE_data.ann$ReferenceID)

q1 <- ggplot(RACE_data.screen, aes(abc_average, fill=TreatmentCode, color=TreatmentCode)) +
  geom_density(alpha=0.25) + 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  facet_wrap(~TreatmentCode) + 
  ggtitle("ABC score distribution at screening") + 
  theme_bw()

q2 <- ggplot(RACE_data.ann, aes(abc_average, fill=TreatmentCode, color=TreatmentCode)) +
  geom_density(alpha=0.25) + 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  facet_wrap(~TreatmentCode) + 
  ggtitle("ABC score distribution at annual follow-up") + 
  theme_bw()

grid.arrange(q1,q2,ncol=1)

#------------------------------------Pairwise t-tests------------------------------------
t.test(
  RACE_data.ann$abc_average[RACE_data.ann$TreatmentCode=="ACTIVESTEP"], 
  RACE_data.screen$abc_average[RACE_data.screen$TreatmentCode=="ACTIVESTEP"], 
  paired = TRUE
)

t.test(
  RACE_data.ann$abc_average[RACE_data.ann$TreatmentCode=="STANDARD"], 
  RACE_data.screen$abc_average[RACE_data.screen$TreatmentCode=="STANDARD"], 
  paired = TRUE
)

#------------------------------------GEE modeling------------------------------------
## Combine complete data from 2 time points:
combined_screen_end <- rbind(RACE_data.screen, RACE_data.ann)
combined_screen_end <- combined_screen_end[order(combined_screen_end$ReferenceID),  ]
length(unique(combined_screen_end$ReferenceID)) #checkpoint: 350 expected

## Convert case-control status to 0,1:
combined_screen_end$TreatmentCode <- gsub("ACTIVESTEP", 1, combined_screen_end$TreatmentCode)
combined_screen_end$TreatmentCode <- gsub("STANDARD", 0, combined_screen_end$TreatmentCode)
combined_screen_end$TreatmentCode <- as.factor(combined_screen_end$TreatmentCode)

## Generate time-order variable:
combined_screen_end$timeorder <- ifelse(combined_screen_end$IntervalName=="Screening", 1, 2)
combined_screen_end$timeorder <- as.integer(combined_screen_end$timeorder)

gee1 <- geeglm(
  abc_average ~ TreatmentCode * timeorder, #note syntax
  data = combined_screen_end, 
  id = ReferenceID, 
  family = poisson("identity"), 
  corstr = "exchangeable"
)

gee1.tab <- as.data.frame(cbind(
  summary(gee1)$coef,
  confint(gee1)
))
colnames(gee1.tab)[colnames(gee1.tab) %in% c("lwr", "upr")] <- c("CI_lower", "CI.upper")
gee1.tab <- round(gee1.tab, 2)
gee1.tab$`Pr(>|W|)`[gee1.tab$`Pr(>|W|)`==0] <- "<2E-16"
grid.arrange(tableGrob(gee1.tab))
