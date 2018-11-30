# Helper functions for XG/RH assessment project
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Date: (ongling)
# Copyright (c) 2018 Y. David Chen
# Notes: 

load_spreadsheet <- function(path, excludeIDs=NULL) {
  #'@description Load each of the 3 worksheets into R
  #'@path path to spreadsheets
  #'@param excludeIDs Accession numbers (rows) to exclude
  data <- read.csv(path, row.names=1, stringsAsFactors=FALSE);
  if(! is.null(excludeIDs)) data <- subset(data, ! rownames(data) %in% excludeIDs); 
  return(data);
}

helper_mcnemar <- function(vec1, vec2) {
  #'@description Helper function for use in iterative Mcnemar's test function
  #'@param vec1,vec2 (Column) Vectors with 0/1
  contTab <- table(vec1, vec2);
  contTab <- contTab[c(2,1), c(2,1)];
  res <- mcnemar.test(contTab, correct=TRUE);
  return(res);
}