# Physician Assessment Concordance Analysis
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Last update: 11/29/2018
# Copyright (c) 2018-19 Y. David Chen
# Notes:

rm(list=ls());
library(pheatmap);
library(psych);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("Helperfunctions_for_RHXG_proj.R");
OUTPUT_XLS <- "~/Downloads/010919_McNemar_tests.xls"; 

analyze_concord <- function(data1, data2, stratMat, adjMet=NULL) {
  #'@description Run independent McNemar's test by patient
  #'@param data1,data2 R matrix or data.frames for data types. Must be matched!
  #'@param stratMat Indicator matrix to stratify (here, it's confidence)
  #'@param adjMet Method for P-value adjustment. Defaults to NULL
  stopifnot(all(rownames(data1) == rownames(data2))); 
  stopifnot(all(rownames(data1) == rownames(stratMat))); 
  resOverall <- resConfid  <- resAmbig <- NULL;
  n <- nrow(data1); 
  
  for(L in colnames(stratMat)) {
    ## Overall:
    res <- helper_twobytwo(data1[ , L], data2[ , L], test=mcnemar.test);
    ck <- psych::cohen.kappa(table(data1[ , L], data2[ , L])); 
    acc <- mean(as.logical(data1[ , L]) == ! as.logical(data2[ , L]) ); #note they are opposite
    resOverall <- rbind(resOverall, c(Proportion=1, Physician=L, Kappa=ck$weighted.kappa, Chisq=res$statistic, Pval=res$p.value, Concordance=acc));
    rm(res);
    
    ##  Confident, non-ambiguous stratum (0):
    obs_confid <- rownames(stratMat)[stratMat[,L] == 0];
    res <- helper_twobytwo(data1[obs_confid, L], data2[obs_confid, L], test=mcnemar.test);
    ck <- psych::cohen.kappa(table(data1[obs_confid, L], data2[obs_confid, L]));
    acc <- mean(data1[obs_confid, L] == data2[obs_confid, L]);
    resConfid <- rbind(resConfid, c(Proportion=length(obs_confid)/n, Physician=L, Kappa=ck$weighted.kappa, Chisq=res$statistic, Pval=res$p.value, Concordance=acc));
    rm(res, ck);
    
    ## Ambiguous stratum (1)
    obs_ambig <- rownames(stratMat)[stratMat[,L] == 1];
    res <- helper_twobytwo(data1[obs_ambig, L], data2[obs_ambig, L], test=mcnemar.test);
    ck <- psych::cohen.kappa(table(data1[obs_ambig, L], data2[obs_ambig, L])); 
    acc <- mean(data1[obs_ambig, L] == data2[obs_ambig, L]);
    resAmbig <- rbind(resAmbig, c(Proportion=length(obs_ambig)/n, Physician=L, Kappa=ck$weighted.kappa, Chisq=res$statistic, Pval=res$p.value, Concordance=acc)); 
  }
  
  ## Final processing & packaging for export:
  resOverall <- as.data.frame(resOverall, stringsAsFactors=FALSE); 
  resConfid <- as.data.frame(resConfid, stringsAsFactors=FALSE);
  resAmbig <- as.data.frame(resAmbig, stringsAsFactors=FALSE);
  if(! is.null(adjMet)) {
    print(paste("Adjusting P-values using the", adjMet, "method..."));
    resOverall$FDR <- p.adjust(resOverall$Pval, method=adjMet);
    resConfid$FDR <- p.adjust(resConfid$Pval, method=adjMet); 
    resAmbig$FDR <- p.adjust(resAmbig$Pval, method=adjMet);
  }
  
  resList <- list(resOverall=resOverall, resConfid=resConfid, resAmbig=resAmbig); 
  return(resList);
}

## Main:
main <- function() {
  ## Load data: 
  xg_data <- load_spreadsheet(paste0(DATA_DIR, "/Column2_Xanthogranuloma.csv"));
  rh_data <- load_spreadsheet(paste0(DATA_DIR, "/Column3_Ret.csv"));
  ambig_indicator <- load_spreadsheet(paste0(DATA_DIR, "/Column4_rater_confidence.csv"));
  # keymap <- load_spreadsheet(paste0(DATA_DIR, "/patient_keymap.csv"));
  
  ## Important checkpoints:
  stopifnot(all(rownames(xg_data) == rownames(rh_data))); 
  stopifnot(all(rownames(xg_data) == rownames(ambig_indicator)));
  
  ## Additional calculations:
  print("Proportion of confident observations:");
  print( 1 - mean(as.matrix(ambig_indicator), na.rm=TRUE) )
  print( prod(dim(ambig_indicator)) - sum(as.matrix(ambig_indicator), na.rm=TRUE) )
  
  ## Inferential statistics:
  resList <- analyze_concord(xg_data, rh_data, ambig_indicator);
  # WriteXLS(resList, ExcelFileName=OUTPUT_XLS); 
  resList
}

main();
