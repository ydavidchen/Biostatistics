# Correlational analysis of clinical variables, part 2
# Script author: David Chen
# Script maintainer: David Chen
# Last update: 04/05/2019
# Copyright (c) 2019 ydavidchen
# Notes:

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions.R");

SCATTER_THEME$axis.title <- element_text(size=8,color="black");
SCATTER_THEME$axis.text <- element_text(size=8,color="black");
SCATTER_THEME$title <- element_text(size=8, color="black", face="bold"); 

sheetWisePipeline <- function(panelCols=5) {
  ## Read data sheet by sheet
  ## Run cor test & plot each sheet, save to plot object + results table
  ## Visualize all at once (consider scaling up canvas)
  plotList <- list();
  for(k in 1:15) {
    dat_k <- load_clin_data(sheet=k);
    plotList[[k]] <- runUnivarLayoutCorTest(dat_k, visualize=TRUE, what2Return="plot.object");
  }
  grid.arrange(grobs=plotList, ncol=panelCols);
}

sheetWisePipeline();
