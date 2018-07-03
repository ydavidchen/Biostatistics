##############################################################################################################
# Pairwise Scatter Plots
# Script author: David Chen
# Script maintainer: David Chen
# Date: 04/02/2018
# Notes:
##############################################################################################################

rm(list=ls())
library(data.table)
library(gdata)
library(GGally)
library(ggplot2)
library(pheatmap)
library(reshape2)

path <- easycsv::choose_dir();
setwd(path); 
recCompl <- read.csv("041918_Cleaned_up_data_for_remodeling.csv", stringsAsFactors=F);
strtDat <- recCompl[ , c("number","strt","Scan_date","N1")];
strtDat$ID <- paste(strtDat$number, strtDat$Scan_date, sep="_");
strtDat$number <- strtDat$Scan_date <- NULL;

## Subset:
strtE <- subset(strtDat, strt=="E");
strtW <- subset(strtDat, strt=="W");
strtN <- subset(strtDat, strt=="N");
strtS <- subset(strtDat, strt=="S");
strtE$strt <- strtW$strt <- strtN$strt <- strtS$strt <- NULL;

strtE <- aggregate(. ~ ID, FUN=mean, data=strtE);
strtW <- aggregate(. ~ ID, FUN=mean, data=strtW);
strtN <- aggregate(. ~ ID, FUN=mean, data=strtN);
strtS <- aggregate(. ~ ID, FUN=mean, data=strtS);

## Further subset:
commonPatients <- Reduce(intersect, list(strtE$ID, strtW$ID, strtN$ID, strtS$ID));
length(commonPatients)

strtE <- subset(strtE, ID %in% commonPatients);
strtW <- subset(strtW, ID %in% commonPatients);
strtN <- subset(strtN, ID %in% commonPatients);
strtS <- subset(strtS, ID %in% commonPatients);

## Merge:
colnames(strtE)[2] <- "E";
colnames(strtW)[2] <- "W";
colnames(strtN)[2] <- "N";
colnames(strtS)[2] <- "S";

strtNWSE <- merge(strtN, strtW, by="ID");
strtNWSE <- merge(strtNWSE, strtS, by="ID");
strtNWSE <- merge(strtNWSE, strtE, by="ID");

## Pairwise scatter plots:
ggpairs(strtNWSE[ , 2:5],
        lower = list(continuous="smooth"),
        diag = list(continuous="barDiag")) +
  theme_classic() +
  theme(axis.text=element_text(size=15,color="black"),
        axis.title=element_text(size=20,color="black"),
        legend.position="top", legend.title=element_blank(),legend.text=element_text(size=20,color="black",face="bold"),
        strip.text=element_text(size=20,colour="black",face="bold"))
