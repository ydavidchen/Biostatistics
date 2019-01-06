# Phase II Clinical Trial Analysis: Overall Survival (OS)
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Copyright (c) 2018 @ydavidchen
# Dates: 06/02/18 (init), 01/06/19 (revision)
# Notes:

rm(list=ls());
if(interactive()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions_folbrite.R");

## Load data:
data <- load_data();

## Initialize storage object for survival results:
osModels <- list();

#--------------------------------------Curve 3: OS based on FLIPI1 (low, int, high) – will be 3 that died--------------------------------------
os_flipi1 <- data;
table(os_flipi1$FLIPI1.low..int..high, useNA="always")
os_flipi1$FLIPI1 <- os_flipi1$FLIPI1.low..int..high;

os_flipi1$FLIPI1 <- gsub("low","Low (N=7)",os_flipi1$FLIPI1); 
os_flipi1$FLIPI1 <- gsub("int","Intermediate (N=17)",os_flipi1$FLIPI1); 
os_flipi1$FLIPI1 <- gsub("high","High (N=15)",os_flipi1$FLIPI1); 

os_flipi1$FLIPI1 <- factor(os_flipi1$FLIPI1, levels=c("Low (N=7)","Intermediate (N=17)","High (N=15)"));
table(os_flipi1$FLIPI1) #check

osModels[["FLIPI1"]] <- survfit(
  Surv(time=os_years, event=os_status) ~ FLIPI1, 
  data = os_flipi1
);

osModels[["FLIPI1"]]

#--------------------------------------Curve 5: OS based on FLIPI2 (low, int, high) – will be 3 that died--------------------------------------
os_flipi2 <- data; 
table(os_flipi2$FLIPI2.low..int..high, useNA="always")
os_flipi2$FLIPI2 <- os_flipi2$FLIPI2.low..int..high;
os_flipi2$FLIPI2 <- gsub("low", "Low (N=6)", os_flipi2$FLIPI2);
os_flipi2$FLIPI2 <- gsub("int", "Intermediate (N=27)", os_flipi2$FLIPI2);
os_flipi2$FLIPI2 <- gsub("high", "High (N=6)", os_flipi2$FLIPI2);

os_flipi2$FLIPI2 <- factor( os_flipi2$FLIPI2, levels=c("Low (N=6)","Intermediate (N=27)","High (N=6)") );
table(os_flipi2$FLIPI2) #check

osModels[["FLIPI2"]] <- survfit(
  Surv(time=os_years, event=os_status) ~ FLIPI2, 
  data = os_flipi2
);

osModels[["FLIPI2"]] 

#--------------------------------------Curve 7: OS based on BCL2 (pos, neg) – will be 3 that died--------------------------------------
os_bcl2 <- data; 
table(os_bcl2$BCL2.positive.at.baseline, useNA="always")
os_bcl2$BCL2[os_bcl2$BCL2.positive.at.baseline=="pos"] <- "Positive (N=12)";
os_bcl2$BCL2[os_bcl2$BCL2.positive.at.baseline=="neg"] <- "Negative (N=27)";
os_bcl2$BCL2 <- factor(os_bcl2$BCL2, levels=c("Positive (N=12)", "Negative (N=27)"));
table(os_bcl2$BCL2) #check

osModels[["BCL2"]] <- survfit(
  Surv(time=os_years, event=os_status) ~ BCL2, 
  data = os_bcl2
);

osModels[["BCL2"]]

#--------------------------------------Data visualization--------------------------------------
for(m in osModels) {
  print( summary(m, time=3) ); #in years
  
  strataName <- gsub( "=.*", "", names(m$strata) );
  strataName <- unique(strataName);
  strataName <- gsub( ".", " ", strataName, fixed=TRUE ); #for legend TITLE
  names(m$strata) <- gsub(".*(=[A-Z])", "\\1", names(m$strata)); #simplify legend TEXT
  names(m$strata) <- substr(names(m$strata), 2, nchar(names(m$strata)));
  
  g <- ggsurvplot(
    m,
    ## Risk Table:
    risk.table.col = "strata",
    risk.table = FALSE,
    tables.theme = theme_cleantable(),
    ## Plot:
    pval = TRUE,
    pval.size = 7,
    conf.int = FALSE,
    conf.int.style = "step",
    censor = TRUE,
    size = 1.5,
    alpha = 0.7,
    ## Labels:
    legend = c(0.75, 0.2), #relative position
    legend.title = strataName,
    xlab = "Time (years)",
    ## Global elements:
    ggtheme = MY_SURV_THEME,
    palette = "Dark2",
    linetype = "strata"
  );
  # print(g)
  ggsave(paste0("~/Downloads/Lansigan_et_al_Figures/OS_",strataName,".png"), dpi=300, width=8, height=8)
}
