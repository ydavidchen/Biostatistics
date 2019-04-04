# Helper functions for clinical correlation analysis
# Script author: David Chen
# Script maintainer: David Chen
# Last update: 04/04/2019
# Copyright (c) 2019 ydavidchen
# Notes:

library(doParallel);registerDoParallel(detectCores() - 1);
library(gdata);
library(ggplot2);
library(gridExtra);

load_clin_data <- function(path="./data/Comparisons for SYNERGY.xlsx", sheet=c(16,17)) {
  #'@description Helper function to load data from multi-tab spreadsheet shared
  #'@param path Path to data
  #'@param sheet Either 16 or 17 (based on data shared)
  if(length(sheet) != 1) stop(paste("You must pick a sheet!, either", sheet));
  data <- read.xls(path, sheet=sheet, stringsAsFactors=FALSE);
  return(data);
}

runCorTests <- function(joinedMat, features, dvName, verbose=TRUE, 
                        visualize=FALSE, plotListName="", panelCols=3) {
  #'@param joinedMat Matrix/dataframe with both independent and dependent variables
  #'@param features Names of independent variables to iterate through
  #'@param dvName Name of dependent variable
  #'@param visualize Optional. Whether multipanel ggplots should be drawn
  #'@param plotListName Optional. Name of object for the plot list, to be returned to global environment
  #'@param panelCols Optional. Number of columns in multi-panel plot layout
  
  resTable <- data.frame(Feature=features, rho=NA, pval=NA, stringsAsFactors=FALSE);
  if(visualize) refinedPlots <- list();
  
  for(feature in features) {
    cTest <- cor.test(joinedMat[,feature], joinedMat[,dvName], method="pearson"); 
    rho <- cTest$estimate;
    pval <- cTest$p.value;
    
    resTable$rho[resTable$Feature==feature] <- rho;
    resTable$pval[resTable$Feature==feature] <- pval;
    
    if(verbose) {
      print(paste("********************", feature, "********************"));
      print(paste("rho:", round(rho, 2) ));
      print(paste("p-value:", signif(pval, 2) )); 
    }
    
    if(visualize) {
      formattedPval <- ifelse(
        pval < 0.01, 
        sprintf("%0.2e", signif(pval, 2)), 
        round(pval, 2)
      ); 
      
      titleLabel <- paste0(
        "rho = ", round(rho, 2), 
        ", p = ", formattedPval 
      );
      
      refinedPlots[[feature]] <- ggplot(joinedMat, aes_string(feature, dvName)) +
        geom_jitter(color="deepskyblue") +
        xlab(gsub(".", " ", feature, fixed=TRUE)) +
        ylab(gsub(".", " ", dvName, fixed=TRUE)) +
        ggtitle(titleLabel) + 
        SCATTER_THEME;
    }
  }
  
  resTable$isSig <- resTable$pval < 0.05;
  
  if(visualize){ 
    grid.arrange(grobs=refinedPlots, ncol=panelCols); #plot instead of return plotList
    if(plotListName != "") assign(plotListName, refinedPlots, envir=.GlobalEnv);
  }
  
  if(verbose) print(resTable);
  return(resTable);
}


SCATTER_THEME <- theme_classic() +
  theme(axis.text=element_text(size=10,color="black"), 
        axis.title=element_text(size=15,color="black"),
        title=element_text(size=12, color="black", face="bold"),
        legend.position="none", legend.title=element_blank(),legend.text=element_text(size=18,color="black"),
        strip.text.x=element_text(size=18,colour="black",face="bold"));

