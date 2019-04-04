# Correlational analysis of clinical variables
# Script author: David Chen
# Script maintainer: David Chen
# Last update: 04/04/2019
# Copyright (c) 2019 ydavidchen
# Notes:

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");

# ---------------------- Manuscript Items 1 - 2 ----------------------
pipeline1617 <- function(y_columns, sheet) {
  #'@description Wrapper to run the common, more comprehensive pipelines
  #'@param y_columns Vector of colums to use as the dependent variable
  #'@param sheet Passed into `load_clin_data` helper
  dat <- load_clin_data(sheet=sheet);
  x_columns <- colnames(dat)[! colnames(dat)  %in% y_columns];
  
  for(dv in y_columns) {
    print(paste("========================", dv, "========================"));
    regr_res <- runCorTests(dat, x_columns, dv, visualize=TRUE);
  }
}

## Item 1:
pipeline1617(c("BMI","Age","Resection.weight"), 16);

## Item 2:
pipeline1617(c("BMI","Age","Resection.Weight"), 17);

# ---------------------- Manuscript Items 3 ----------------------
pipelineSelected <- function(y_columns, x_columns, plotListName, panelCols) {
  #'@description Wrapper to run more customized correlation analyses
  dat <- load_clin_data(sheet=17);
  
  for(dv in y_columns) {
    print(paste("========================", dv, "========================"));
    regr_res <- runCorTests(dat, x_columns, dv, visualize=TRUE, plotListName=plotListName, panelCols=panelCols);
  }
}

pipelineSelected(
  y_columns = c("Satisfaction.POD3"), 
  x_columns = c( "NSAID.dosage.POD3","Pain.scale.POD3","Acetaminophen.dosage.POD3"),
  plotListName = "plotColGs",
  panelCols = 3
);

pipelineSelected(
  y_columns = c("Pain.scale.POD3"), 
  x_columns = c("ME.used.POD3"), 
  plotListName = "plotColE", 
  panelCols = 1
);

combinedList <- plotColGs;
combinedList[["ME.prescribed.POD3"]] <- plotColE$ME.used.POD3;
grid.arrange(grobs=combinedList, ncol=2);
