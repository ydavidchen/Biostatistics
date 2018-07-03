##############################################################################################################
# Line Plot over Time (days)
# Script author: David Chen
# Script maintainer: David Chen
# Date: 04/23/2018
# Notes:
##############################################################################################################

rm(list=ls())
library(data.table)
library(gdata)
library(ggplot2)
library(lubridate)
library(reshape2)

## Load data:
easycsv::choose_dir();
recCompl <- read.csv("041918_Cleaned_up_data_for_remodeling.csv", stringsAsFactors=F);
strtDat <- recCompl[ , c("number","strt","Scan_date","N1")];
strtDat$Day_abs <- as.numeric(as.Date(strtDat$Scan_date));

## Compute relative date:
## The first record on file is defined as day 0
strtDat$Day_relative <- NA;
for(p in unique(strtDat$number)) {
  charInd <- which(strtDat$number == p); #indexing based on character (patient)
  origin_p <- as.numeric(as.Date(min(strtDat$Scan_date[charInd]))); #1st day for a given patient
  strtDat$Day_relative[charInd] <- strtDat$Day_abs[charInd] - origin_p; #diff in days
}
strtDat$ID <- paste(strtDat$number, strtDat$Day, sep="_");

pKeep <- as.data.frame(table(strtDat$number));
pKeep <- as.character(pKeep$Var1[pKeep$Freq > 4])

## Data visualization
ggplot(subset(strtDat, number %in% pKeep), aes(x=Day_relative, y=N1, color=ID)) +
  geom_line() +
  labs(x="Number of days (relative to patient's first record)", y="Penetration distance (mm)") +
  theme_bw() +
  theme(axis.text=element_text(size=16,color="black"), axis.title=element_text(size=16,color="black"),
        legend.position="none", legend.title=element_blank(),legend.text=element_text(size=14,color="black",face="bold"),
        strip.text.x=element_text(size=12,colour="black",face="bold")) +
  facet_wrap(~strt)
