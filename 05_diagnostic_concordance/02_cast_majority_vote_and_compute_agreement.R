# Assign Diagnosis by Majorty Vote 
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Last update: 01/25/2019
# Copyright (c) 2018-19 Y. David Chen
# Notes:

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("Helperfunctions_for_RHXG_proj.R");
OUTPUT_NAME <- "~/Downloads/012519_sum_data_and_stat_results.xls";

test_confidence_against_majority_vote <- function(newDf, outcome, observers) {
  dataVec1 <- newDf[ , outcome];
  resDf <- data.frame(Observer=observers, oddsRatio=NA, ci.lb=NA, ci.ub=NA, pval=NA); #init
  for(j in observers){
    dataVec2 <- newDf[ , j];
    test_res <- helper_twobytwo(dataVec1, dataVec2, testFUN=fisher.test);
    resDf$oddsRatio[resDf$Observer == j] <- test_res$estimate;
    resDf$pval[resDf$Observer == j] <- test_res$p.value;
    resDf[resDf$Observer == j, c("ci.lb","ci.ub")] <- test_res$conf.int;
  }
  resDf$isSignif <- resDf$pval <= 0.05;
  return(resDf);
}

calcMode <- function(vec) {
  #'@description Helper function to calculate mode of a numeric array. Passed into `apply`
  #'@param vec Numerical vector
  uniqv <- unique(vec); 
  res <- uniqv[which.max(tabulate(match(vec, uniqv)))];
  return(res); 
}

agreement_fun <- function(vec) {
  #'@description Helper function to calculate %agreement aka concordance. Paseed into `apply`
  #'@param vec Numerical vector
  mu_1 <- mean(vec == 1);
  mu_0 <- mean(vec == 0);
  res <- max(mu_1, mu_0); 
  return(res);
}

main <- function() {
  ## Load data: 
  xg_data <- load_spreadsheet(paste0(DATA_DIR, "/Column2_Xanthogranuloma.csv"));
  rh_data <- load_spreadsheet(paste0(DATA_DIR, "/Column3_Ret.csv"));
  
  ambig_indicator <- load_spreadsheet(paste0(DATA_DIR, "/Column4_rater_confidence.csv"));
  conf_indicator <- ifelse(ambig_indicator==1, 0, 1);
  stopifnot(all(ambig_indicator+conf_indicator==1)); #checkpoint
  
  # ---------------------------------- Majority Diagnosis vs. Confidence ---------------------------------- 
  ## Cast majority vote:
  xg_data$isXG <- apply(xg_data, 1, calcMode); 
  rh_data$isRH <- apply(rh_data, 1, calcMode); 
  
  ## Statistical tests:
  newDf <- cbind(
    xg_data[,"isXG", drop=FALSE], 
    rh_data[,"isRH", drop=FALSE]
  );
  newDf <- merge(newDf, conf_indicator, by="row.names");
  stopifnot(all(newDf$isXG + newDf$isRH == 1)); #check for mutual exclusivity
  statRes <- test_confidence_against_majority_vote(newDf, "isXG", colnames(conf_indicator));
  
  ## Data visualization:
  colnames(xg_data)[colnames(xg_data)=="isXG"] <- "Diagnosis";
  colnames(rh_data)[colnames(rh_data)=="isRH"] <- "Diagnosis";
  
  custom_heatmap(xg_data, main="XG", gaps_col=5);
  custom_heatmap(rh_data, main="RH", gaps_col=5);
  custom_heatmap(ambig_indicator, main="Self-reported Ambiguous \n (1=ambiguous, 0=confident)");
  
  # ---------------------------------- Descriptive Report of %Agreement ---------------------------------- 
  ## Compute probability (similar idea as majority vote):
  colNames <- LETTERS[1:5]; 
  xg_data$Agreement <- apply(xg_data[,colNames], 1, FUN=agreement_fun);
  rh_data$Agreement <- apply(rh_data[,colNames], 1, FUN=agreement_fun);
  WriteXLS::WriteXLS(list(res=statRes,XG=xg_data,RH=rh_data), ExcelFileName=OUTPUT_NAME, row.names=TRUE);
}

main();
