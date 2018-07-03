##############################################################################################################
# Regression Modeling
# Script author: David Chen
# Script maintainer: David Chen
# Date: 06/07/2018
# Notes:
##############################################################################################################

rm(list=ls())
library(lme4)
library(lmerTest)
library(MuMIn)
library(doParallel); registerDoParallel(detectCores() - 1)

## Load data:
easycsv::choose_dir()
recCompl <- read.csv("FilterPerf_legs060118_cleaned.csv", stringsAsFactors=FALSE); 

## Time series regression:
fitRE <- lmer(
  N1 ~ years*sex + Age_at_first_proc + sex + dIVC + Thrombosis + Anticoagulation + Cancer + (1 | number / strt), 
  data = recCompl
);
s <- summary(fitRE); #P-value requires lmerTest
s

## Polynomial time series:
fitPolyRE <- lmer(
  N1 ~ poly(years, 2)*sex + Age_at_first_proc + sex + dIVC + Thrombosis + Anticoagulation + Cancer + (1 | number / strt), 
  data = recCompl
);
sPoly <- summary(fitPolyRE); #P-value requires lmerTest
sPoly

## Logarithmic time series:
fitLogRE <- lmer(
  N1 ~ log(years + 1/365.25)*sex + Age_at_first_proc + sex + dIVC + Thrombosis + Anticoagulation + Cancer + (1 | number / strt), 
  data = recCompl
);
sLog <- summary(fitLogRE); #P-value requires lmerTest
sLog

## Likelihood Ratio Test:
anova(fitRE, fitPolyRE, test="Chisq")
anova(fitRE, fitLogRE, test="Chisq")

## Summary of pseudo R-squared
pseudoR_sq <- rbind(
  r.squaredGLMM(fitRE),
  r.squaredGLMM(fitPolyRE),
  r.squaredGLMM(fitLogRE)
);
rownames(pseudoR_sq) <- c("Linear", "Polynomial", "Logarithmic");
pseudoR_sq 
