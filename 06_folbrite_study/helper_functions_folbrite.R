# Phase II Clinical Trial Analysis: helper functions
# Script author: David Chen
# Script maintainer: David Chen
# Date: (ongoing)
# Notes:

library(ggplot2);
library(doParallel); registerDoParallel(detectCores() - 1);

load_data <- function(path="~/Dropbox (Christensen Lab)/Christensen Lab - 2018/StatsConsulting_2018/Survival_curve_proj/Sheet1_Updated.csv") {
  df <- read.csv(path, stringsAsFactors=FALSE);
  df[df == ""] <- NA; 
  stopifnot( ! any(duplicated(df$Patient.ID)) ); 
  return(df);
}

sentenceCase <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1));
  return(x);
}

## Graphic themes:
mySurvTheme <- theme_classic2() +
    theme(axis.text=element_text(size=20,color="black"), 
          axis.title=element_text(size=20,color="black"),
          legend.title=element_text(size=16,color="black",face="bold"), legend.text=element_text(size=16,color="black",face="bold"), 
          legend.box.background=element_rect(colour="black",size=2),
          text=element_text(size=20),
          strip.text.x=element_text(size=15,colour="black",face="bold"));
  
