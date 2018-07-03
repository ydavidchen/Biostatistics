##############################################################################################################
# Data Cleaning & Exploration
# Script author: David Chen
# Script maintainer: David Chen
# Date: 03/15/2018; 04/01/2018
# Notes:
##############################################################################################################

rm(list=ls())
library(plyr); library(dplyr)
library(gdata)
library(ggplot2)
library(lme4); library(lmerTest)
library(reshape2)
library(tableone)

## Load cleaned, resaved data:
easycsv::choose_dir()
patientRec <- read.csv("FilterPerf_legs030118_cleaned.csv", stringsAsFactors=F)
patientRec[patientRec == ""] <- NA;

#------------------------------------Table One------------------------------------
## Select unique records for table one and dropping:
recUniq <- patientRec[! duplicated(patientRec$number), ];

myFactorVars <- c("Cancer", "Thrombosis", "Anticoagulation", "sex", "VTE", "Femoral_access");
myVars <- c("Age", "initial_diameter", myFactorVars);
myTableOne <- CreateTableOne(
  data = recUniq, 
  vars = myVars, 
  factorVars = myFactorVars,
  smd = FALSE, 
  test = FALSE,
  includeNA = TRUE
);
y <- print(myTableOne, showAllLevels=TRUE); 
# write.csv(y, file="~/Downloads/031518_Patient_records_table1.csv", row.names=T, quote=F)

#------------------------------------Random effect modeling using mean of 4 filter legs------------------------------------
recCompl <- patientRec; #initialize
recCompl <- recCompl[ , c("number","Scan_date","strt","N1")];
recCompl <- subset(recCompl, ! is.na(Scan_date));
rownames(recCompl) <- NULL;

## Add temporal order:
recCompl$Scan_date <- as.Date(factor(recCompl$Scan_date));
recCompl$days <- NA;
for(n in unique(recCompl$number)) {
  day1 <- recCompl$Scan_date[recCompl$number == n][1];
  recCompl$days[recCompl$number==n] <- recCompl$Scan_date[recCompl$number == n] - day1;
}
hist(recCompl$days)
if(anyNA(recCompl$days)) {
  print("Dropping rows with missing dates")
  recCompl <- subset(recCompl, ! is.na(days) ); 
} else {
  print("No missing value in dates.")
}
recCompl$years <- recCompl$days / 365;

## Merge in covariates:
recCompl <- merge(
  recCompl,
  recUniq[ , c("number","Age","sex","dIVC","Cancer","Thrombosis","Anticoagulation")], 
  by = "number", 
  all.x = TRUE
);
# write.csv(recCompl, file="~/Dropbox (Christensen Lab)/Christensen Lab - 2018/QBS123_2018/Cardiology_Proj/041918_Cleaned_up_data_for_remodeling.csv", row.names=F, quote=F)

## Time series regression accounting for autocorrelation (random effect):
fitRE <- lmer(
  N1 ~ years + (1 | number / strt) + Age + sex + dIVC + Thrombosis + Anticoagulation + Cancer, 
  data = recCompl
);
s <- summary(fitRE) #P-value requires lmerTest
s

#-------------------------------------------------Data visualization-------------------------------------------------
## Customized forest plots:
coefMat <- as.data.frame( s$coefficients );
coefMat$lower <- coefMat$Estimate - coefMat$`Std. Error`;
coefMat$upper <- coefMat$Estimate + coefMat$`Std. Error`;
coefMat$Category <- rownames(coefMat);
coefMat <- subset(coefMat, Category != "(Intercept)");
coefMat$Category <- gsub("years", "Time (years)", coefMat$Category);
coefMat$Category <- gsub("days", "Time (days)", coefMat$Category);
coefMat$Category <- gsub("Age", "Age (years)", coefMat$Category);
coefMat$Category <- gsub("sexm", "Sex (male)", coefMat$Category);
coefMat$Category <- gsub("TRUE", "", coefMat$Category);
coefMat$Category <- factor(coefMat$Category, levels=rev(coefMat$Category)); #use as-is levels 
coefMat$P <- signif(coefMat$`Pr(>|t|)`, 3);
ggplot(coefMat, aes(x=Category, y=Estimate, ymin=lower, ymax=upper)) +
  geom_pointrange(size=1) +
  geom_hline(yintercept=0, size=0.3, linetype="dashed") +
  coord_flip() + #order: left to right becomes bottom to up
  scale_y_continuous(limits=c(-0.7,0.5)) + #consistent scale for all
  theme_classic() +
  theme(axis.line=element_line(color="black"), axis.ticks=element_line(color="black"),
        axis.text=element_text(size=20, color="black"), title=element_text(size=10, color="black",face="italic"),
        axis.title.x=element_text(size=20,face="bold"), axis.title.y=element_blank()) +
  labs(x="Fixed-effect coefficients", y="Coefficient estimate") +
  ggtitle("Distance ~ Time (years) + (1 | patient | strt) + Age + Sex + Baseline (dIVC) + Thrombosis + Anticoagulation + Cancer") + 
  # annotate("text", 6.15, coefMat$Estimate[1], label=paste("P =", coefMat$P[1]),size=7) +
  annotate("text", 7.15, coefMat$Estimate[1], label="P < 2e-16",size=7) +
  annotate("text", 6.15, coefMat$Estimate[2], label=paste("P =", coefMat$P[2]),size=7) +
  annotate("text", 5.15, coefMat$Estimate[3], label=paste("P =", coefMat$P[3]),size=7) +
  annotate("text", 4.15, coefMat$Estimate[4], label=paste("P =", coefMat$P[4]),size=7) +
  annotate("text", 3.15, coefMat$Estimate[5], label=paste("P =", coefMat$P[5]),size=7) + 
  annotate("text", 2.15, coefMat$Estimate[6], label=paste("P =", coefMat$P[6]),size=7) +
  annotate("text", 1.15, coefMat$Estimate[7], label=paste("P =", coefMat$P[7]),size=7)

## Scatter plot of raw data points:
ggplot(recCompl, aes(x=days, y=Mean)) +
  geom_point(aes(color=sex)) +
  geom_smooth(method="lm", se=FALSE, color="black") + 
  theme_classic() +
  theme(axis.text.x=element_text(size=15,color="black"), axis.title.x=element_text(size=15,color="black"),
        axis.text.y=element_text(size=15,color="black"), axis.title.y=element_text(size=15,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"), title=element_text(size=10, color="black",face="italic"),
        legend.position="top", legend.title=element_text(size=15,color="black"), legend.text=element_text(size=15,color="black")) +
  ggtitle("Model: Mean of 4 filters ~ Time (days) + Age + Sex + Thrombosis + Anticoagulation + Cancer + Subject (random effect)") + 
  annotate(geom="text",x=4000,y=5.5,label="P = 8.88e-16 ***",size=6) +
  annotate(geom="text",x=4000,y=5,label="Coefficient = 3.8e-04",size=6)

