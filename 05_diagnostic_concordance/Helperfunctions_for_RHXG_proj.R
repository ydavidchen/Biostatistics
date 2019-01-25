# Helper functions for XG/RH assessment project
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Date: (ongling)
# Copyright (c) 2018 Y. David Chen
# Notes: 

library(doParallel); registerDoParallel(detectCores() - 1);
library(WriteXLS);
DATA_DIR <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/StatsConsulting_2018/RHXG_study/rhxg_data/";

load_spreadsheet <- function(path, excludeIDs="SD-16-26909", excludePhys="F") {
  #'@description Load each of the 3 worksheets into R
  #'@path path to spreadsheets
  #'@param excludeIDs Accession numbers (rows) to exclude
  #'@param excludePhys Observer to exclude
  data <- read.csv(path, row.names=1, stringsAsFactors=FALSE);
  if(! is.null(excludeIDs)) data <- subset(data, ! rownames(data) %in% excludeIDs); 
  if(! is.null(excludePhys)) data <- data[ , ! colnames(data) %in% excludePhys]; 
  return(data);
}

custom_heatmap <- function(data, ...) {
  #'@description Customized pheatmap for consistency
  #'@param data Numerical data matrix to plot
  #'@param ... Additional parameters passed to pheatmap
  require(pheatmap);
  ph <- pheatmap(
    data,  
    color = c("lightgray","black"),
    border_color = NA,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_col = 13, 
    legend_breaks = c(0,1),
    ...
  ); 
  return(ph);
}

helper_twobytwo <- function(vec1, vec2, testFUN) {
  #'@description Helper function for use in iterative two-by-two discrete tests
  #'@param vec1,vec2 (Column) Vectors with 0/1
  #'@param testFUN Either `fisher.test` or `mcnemar.test`
  contTab <- table(vec1, vec2);
  contTab <- contTab[c(2,1), c(2,1)];
  print(contTab);
  res <- testFUN(contTab);
  return(res);
}

