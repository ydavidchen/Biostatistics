# Response to Reviewer Question about ad hoc Power Calculation
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Date: 11/28/2018
# Notes:
# -- Be sure to state the assumptions!

rm(list=ls());
library(pwr); 

SIG_THRESH <- 0.05; 
PWR_THRESH <- 0.80; 

N_A <- 236;
MU_A <- 9.56;
SD_A <- 1.30;

N_B <- 404; 
MU_B <- 9.59;
SD_B <- 0.95; 

pooledSD <- (SD_A + SD_B) / 2;
pooledSD

ptt <- pwr.t.test(
  n = (N_A+N_B) / 2,
  d = NULL,
  sig.level = SIG_THRESH, 
  power = PWR_THRESH, 
  type = "one.sample", 
  alternative = "two.sided"
);
ptt

sol <- ptt$d * pooledSD;
paste("Expected FC:", sol)
paste("Mean difference detectable", MU_A*sol)
