# Physician Assessment Concordance Analysis
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Last update: 11/29/2018
# Copyright (c) 2018 Y. David Chen
# Notes: 
# -- Run this script in RStudio

rm(list=ls());
library(psych);
library(pheatmap);
fixInNamespace("draw_colnames","pheatmap");  #vjust = 1, hjust = 0.5, rot = 0
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("Helperfunctions_for_RHXG_proj.R");
DATA_DIR <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/StatsConsulting_2018/RHXG_study/rhxg_data/";
EXCLUDED <- c("SD-16-26909");

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
    ...
  ); 
  return(ph);
}

analyze_concord <- function(data1, data2, stratMat, adjMet="fdr", sigThresh=0.05) {
  #'@description Run independent McNemar's test by patient
  #'@param data1,data2 R matrix or data.frames for data types. Must be matched!
  #'@param stratMat Indicator matrix to stratify (here, it's confidence)
  stopifnot(all(rownames(data1) == rownames(data2))); 
  stopifnot(all(rownames(data1) == rownames(stratMat))); 
  resOverall <- resConfid  <- resAmbig <- NULL;
  n <- nrow(data1); 
  
  for(L in colnames(stratMat)) {
    ## Overall:
    res <- helper_mcnemar(data1[ , L], data2[ , L]);
    ck <- psych::cohen.kappa(table(data1[ , L], data2[ , L])); 
    resOverall <- rbind(resOverall, c(Proportion=1, Physician=L, Kappa=ck$weighted.kappa, Chisq=res$statistic, Pval=res$p.value));
    rm(res);
    
    ## Confident stratum:
    obs_confid <- rownames(self_conf)[stratMat[,L] == 1];
    res <- helper_mcnemar(data1[obs_confid, L], data2[obs_confid, L]);
    ck <- psych::cohen.kappa(table(data1[obs_confid, L], data2[obs_confid, L]));
    resConfid <- rbind(resConfid, c(Proportion=length(obs_confid)/n, Physician=L, Kappa=ck$weighted.kappa, Chisq=res$statistic, Pval=res$p.value));
    rm(res, ck);
    
    ## Ambiguous stratum:
    obs_ambig <- rownames(self_conf)[stratMat[,L] == 0];
    res <- helper_mcnemar(data1[obs_ambig, L], data2[obs_ambig, L]);
    ck <- psych::cohen.kappa(table(data1[obs_ambig, L], data2[obs_ambig, L])); 
    resAmbig <- rbind(resAmbig, c(Proportion=length(obs_ambig)/n, Physician=L, Kappa=ck$weighted.kappa, Chisq=res$statistic, Pval=res$p.value)); 
  }
  
  ## Final processing & packaging for export:
  resOverall <- as.data.frame(resOverall, stringsAsFactors=FALSE); 
  resOverall$FDR <- p.adjust(resOverall$Pval, method=adjMet);
  resOverall$isSignif <- resOverall$FDR <= sigThresh; 
  
  resConfid <- as.data.frame(resConfid, stringsAsFactors=FALSE);
  resConfid$FDR <- p.adjust(resConfid$Pval, method=adjMet); 
  resConfid$isSignif <- resConfid$FDR <= sigThresh;
  
  resAmbig <- as.data.frame(resAmbig, stringsAsFactors=FALSE);
  resAmbig$FDR <- p.adjust(resAmbig$Pval, method=adjMet);
  resAmbig$isSignif <- resAmbig$FDR <= sigThresh;
  
  resList <- list(resOverall=resOverall, resConfid=resConfid, resAmbig=resAmbig); 
  return(resList);
}

## Main:
## Load data: 
xg_data <- load_spreadsheet(paste0(DATA_DIR, "/Column2_Xanthogranuloma.csv"), excludeIDs=EXCLUDED);
rh_data <- load_spreadsheet(paste0(DATA_DIR, "/Column3_Ret.csv"), excludeIDs=EXCLUDED);
self_conf <- load_spreadsheet(paste0(DATA_DIR, "/Column4_rater_confidence.csv"), excludeIDs=EXCLUDED);
# keymap <- load_spreadsheet(paste0(DATA_DIR, "/patient_keymap.csv"));

## Important checkpoints:
stopifnot(all(rownames(xg_data) == rownames(rh_data))); 
stopifnot(all(rownames(xg_data) == rownames(self_conf)));

## Additional calculations:
print( 1 - mean(as.matrix(self_conf), na.rm=TRUE) )
print( prod(dim(self_conf)) - sum(as.matrix(self_conf), na.rm=TRUE) )

## Data visualization:
custom_heatmap(xg_data, main="XG");
custom_heatmap(rh_data, main="RH");
custom_heatmap(self_conf, main="Self-reported Confidence");

## Inferential statistics:
## Consider stratified analysis
resList <- analyze_concord(xg_data, rh_data, self_conf);
resList
WriteXLS::WriteXLS(resList, ExcelFileName="~/Downloads/113018_McNemar_tests.xls"); 
