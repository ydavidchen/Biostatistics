# Phase II Clinical Trial Analysis: Progression Free Survival
# Script author: Y. David Chen
# Script maintainer: Y. David Chen
# Date: 06/01/18
# Notes:

rm(list=ls()); 
library(ggsci);
library(stringr);
library(survival);
library(survminer);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("helper_functions_folbrite.R");

## Load data:
data <- load_data();
data <- subset(data, ! is.na(PFS.Days..n.36.)); 
data$years <- data$PFS.Days..n.36. / 365.25;

## Initialize storage object:
pfsModels <- list();

#--------------------------------------Curve 1: PFS based on post BR response (CR/CRu vs. PR)--------------------------------------
pfs_fu <- data;
table(pfs_fu$PostBR, useNA="always")
pfs_fu$Post.BR[pfs_fu$PostBR %in% c("CRu", "CR")] <- "CR/CRu (N=20)";
pfs_fu$Post.BR[pfs_fu$PostBR %in% c("PR", "SD")] <- "PR (N=16)"; 
table(pfs_fu$Post.BR) #check
pfs_fu$Post.BR <- as.factor(pfs_fu$Post.BR); 

pfsModels[["postBR"]] <- survfit(
  Surv(time=years, event=Event.Free.vs..Progression..n.39.) ~ Post.BR, 
  data = pfs_fu
);
pfsModels[["postBR"]]

## Question from Darcie/Erick:
summary(pfsModels[["postBR"]], times=2) #24months
summary(pfsModels[["postBR"]], times=2.5) #30months

#--------------------------------------Curve 2: PFS based on FLIPI1 (low, int, high)--------------------------------------
pfs_flipi1 <- data;
table(pfs_flipi1$FLIPI1.low..int..high)

pfs_flipi1$FLIPI1 <- pfs_flipi1$FLIPI1.low..int..high;
pfs_flipi1$FLIPI1 <- gsub("low", "Low (N=6)", pfs_flipi1$FLIPI1); 
pfs_flipi1$FLIPI1 <- gsub("int", "Intermediate (N=16)", pfs_flipi1$FLIPI1); 
pfs_flipi1$FLIPI1 <- gsub("high", "High (N=14)", pfs_flipi1$FLIPI1); 

pfs_flipi1$FLIPI1 <- factor(pfs_flipi1$FLIPI1, levels=c("Low (N=6)","Intermediate (N=16)","High (N=14)")); 
table(pfs_flipi1$FLIPI1)

pfsModels[["FLIPI1"]] <- survfit(
  Surv(time=years, event=Event.Free.vs..Progression..n.39.) ~ FLIPI1, 
  data = pfs_flipi1
);
pfsModels[["FLIPI1"]]

#--------------------------------------Curve 4: PFS based on FLIPI2 (low, int, high)--------------------------------------
pfs_flipi2 <- data;
table(pfs_flipi1$FLIPI2.low..int..high)

pfs_flipi2$FLIPI2 <- pfs_flipi1$FLIPI2.low..int..high;
pfs_flipi2$FLIPI2 <- gsub("low", "Low (N=5)", pfs_flipi2$FLIPI2);
pfs_flipi2$FLIPI2 <- gsub("int", "Intermediate (N=25)", pfs_flipi2$FLIPI2);
pfs_flipi2$FLIPI2 <- gsub("high", "High (N=6)", pfs_flipi2$FLIPI2);

pfs_flipi2$FLIPI2 <- factor(pfs_flipi2$FLIPI2, c("Low (N=5)","Intermediate (N=25)","High (N=6)"));
table(pfs_flipi2$FLIPI2)

pfsModels[["FLIPI2"]] <- survfit(
  Surv(time=years, event=Event.Free.vs..Progression..n.39.) ~ FLIPI2, 
  data = pfs_flipi2
);
pfsModels[["FLIPI2"]]

#--------------------------------------Curve 6: PFS based on BCL2 (pos, neg)--------------------------------------
pfs_bcl2 <- data;
table(pfs_bcl2$BCL2.positive.at.baseline)

pfs_bcl2$BCL2 <- pfs_bcl2$BCL2.positive.at.baseline;
pfs_bcl2$BCL2 <- gsub("neg", "Negative (N=26)", pfs_bcl2$BCL2);
pfs_bcl2$BCL2 <- gsub("pos", "Positive (N=10)", pfs_bcl2$BCL2);

pfs_bcl2$BCL2 <- factor(pfs_bcl2$BCL2, levels=c("Positive (N=10)","Negative (N=26)"));

pfsModels[["BCL2"]] <- survfit(
  Surv(time=years, event=Event.Free.vs..Progression..n.39.) ~ BCL2,
  data = pfs_bcl2
);
pfsModels[["BCL2"]] 

#--------------------------------------Curve 8: PFS based on BCL2 pos/neg/FLIPI1 (low/int/high) – 6 curves--------------------------------------
pfs_bcl2xflip1 <- data;
pfs_bcl2xflip1$BCL2.and.FLIPI1 <- paste(pfs_bcl2xflip1$BCL2.positive.at.baseline, pfs_bcl2xflip1$FLIPI1.low..int..high, sep=" & ");
table(pfs_bcl2xflip1$BCL2.and.FLIPI1)

pfs_bcl2xflip1$BCL2.and.FLIPI1 <- gsub("neg & high", "BCL2-, FLIPI1 high (N=10)", pfs_bcl2xflip1$BCL2.and.FLIPI1); 
pfs_bcl2xflip1$BCL2.and.FLIPI1 <- gsub("neg & int",  "BCL2-, FLIPI1 int. (N=11)", pfs_bcl2xflip1$BCL2.and.FLIPI1); 
pfs_bcl2xflip1$BCL2.and.FLIPI1 <- gsub("neg & low",  "BCL2-, FLIPI1 low  (N=5)", pfs_bcl2xflip1$BCL2.and.FLIPI1); 
pfs_bcl2xflip1$BCL2.and.FLIPI1 <- gsub("pos & high", "BCL2+, FLIPI1 high (N=4)", pfs_bcl2xflip1$BCL2.and.FLIPI1); 
pfs_bcl2xflip1$BCL2.and.FLIPI1 <- gsub("pos & int",  "BCL2+, FLIPI1 int. (N=5)", pfs_bcl2xflip1$BCL2.and.FLIPI1); 
pfs_bcl2xflip1$BCL2.and.FLIPI1 <- gsub("pos & low",  "BCL2+, FLIPI1 low  (N=1)", pfs_bcl2xflip1$BCL2.and.FLIPI1); 

table(pfs_bcl2xflip1$BCL2.and.FLIPI1)

pfsModels[["BCL2xFLIPI1"]] <- survfit(
  Surv(time=years, event=Event.Free.vs..Progression..n.39.) ~ BCL2.and.FLIPI1, 
  data = pfs_bcl2xflip1
);
pfsModels[["BCL2xFLIPI1"]]

#--------------------------------------Curve 9: PFS based on BCL2 (pos/neg)/FLIPI2 (low/int/high) – 6 curves--------------------------------------
pfs_bcl2xflip2 <- data;
pfs_bcl2xflip2$BCL2.and.FLIPI2 <- paste(pfs_bcl2xflip2$BCL2.positive.at.baseline, pfs_bcl2xflip2$FLIPI2.low..int..high, sep=" & "); 
table(pfs_bcl2xflip2$BCL2.and.FLIPI2, useNA="ifany")

pfs_bcl2xflip2$BCL2.and.FLIPI2 <- gsub("neg & high", "BCL2-, FLIPI2 high (N=3)", pfs_bcl2xflip2$BCL2.and.FLIPI2);
pfs_bcl2xflip2$BCL2.and.FLIPI2 <- gsub("neg & int",  "BCL2-, FLIPI2 int. (N=18)", pfs_bcl2xflip2$BCL2.and.FLIPI2);
pfs_bcl2xflip2$BCL2.and.FLIPI2 <- gsub("neg & low",  "BCL2-, FLIPI2 low  (N=5)", pfs_bcl2xflip2$BCL2.and.FLIPI2);
pfs_bcl2xflip2$BCL2.and.FLIPI2 <- gsub("pos & high", "BCL2+, FLIPI2 high (N=3)", pfs_bcl2xflip2$BCL2.and.FLIPI2);
pfs_bcl2xflip2$BCL2.and.FLIPI2 <- gsub("pos & int",  "BCL2+, FLIPI2 int. (N=7)", pfs_bcl2xflip2$BCL2.and.FLIPI2);

table(pfs_bcl2xflip2$BCL2.and.FLIPI2 )

pfsModels[["BCL2&FLIPI2"]] <- survfit(
  Surv(time=years, event=Event.Free.vs..Progression..n.39.) ~ BCL2.and.FLIPI2, 
  data = pfs_bcl2xflip2
);
pfsModels[["BCL2&FLIPI2"]]

#--------------------------------------Curve 11: Post-RIT--------------------------------------
pfs_rit <- subset(data, PostRIT != "PD");
table(pfs_rit$PostRIT, useNA="ifany")
pfs_rit$Post.RIT[pfs_rit$PostRIT %in% c("CR","CRu")] <- "CR/CRu (N=30)";
pfs_rit$Post.RIT[pfs_rit$PostRIT == "PR"] <- "PR (N=5)";
table(pfs_rit$Post.RIT, useNA="ifany")

pfsModels[["postRIT"]] <- survfit(
  Surv(time=years, event=Event.Free.vs..Progression..n.39.) ~ Post.RIT, 
  data = pfs_rit
);

#--------------------------------------Data visualization--------------------------------------
for(m in pfsModels) {
  # print( summary(m, time=3) ) #in years
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
    legend = c(0.7, 0.2), #rel. position
    legend.title = strataName,
    xlab = "Time (years)",
    ## Global elements:
    ggtheme = mySurvTheme,
    palette = "Dark2",
    linetype = "strata"
  );
  # print(g)
  ggsave(paste0("~/Downloads/Lansigan_et_al_Figures/PFS_",strataName,".png"), dpi=300, width=8, height=8)
}

rm(m, g);
