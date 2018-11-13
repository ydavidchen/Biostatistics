# Response to Reviewer Comment about Multiple Hypothesis Testing (R. Barth group manuscript)
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Date: 11/13/2018
# Note:
# -- Each table received was resaved as a separate CSV file in Excel

rm(list=ls());
library(data.table);
library(WriteXLS);

ADJ_MET <- "fdr"; 
SIG_THRESH <- 0.05; 
TAB_NAMES <- c("table1", "table2", "table3"); 
TABLE_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/StatsConsulting_2018/Surgeon_manuscript/";
outPath <- paste0(TABLE_PATH, "/outputs");

setwd(TABLE_PATH);

## Table 1:
## This table has 2 analyses
table1 <- fread("table1_resaved.csv", stringsAsFactors=FALSE);

table1$tTestFDR <- p.adjust(table1$`T- Test`, method=ADJ_MET); 
table1$isTtestSignif <- table1$tTestFDR < SIG_THRESH;

table1$ChiSqFDR <- p.adjust(as.numeric(gsub("< ", "", table1$`Chi Squared`)), method=ADJ_MET); #initialize
table1$isChisqSignif <- table1$ChiSqFDR < SIG_THRESH;

## Table 2:
table2 <- fread("table2_resaved.csv", stringsAsFactors=FALSE);
nSkip <- sum(is.na(table2$`T-Test`)); #1
table2$tTestFDR <- p.adjust(table2$`T-Test`, method=ADJ_MET, n=nrow(table2)-nSkip);
table2$isTtestSignif <- table2$tTestFDR < SIG_THRESH;

## Table 3: 
table3 <- fread("table3_resaved.csv", stringsAsFactors=FALSE);
table3$tTestFDR <- p.adjust(as.numeric(gsub("< ","",table3$`T- Test`)), method=ADJ_MET); 
table3$isTTestSignif <- table3$tTestFDR < SIG_THRESH;

## Export:
WriteXLS(
  x = TAB_NAMES, 
  ExcelFileName = paste0(outPath, "/Tables_with_FDR.xls"),
  SheetNames = TAB_NAMES
);
