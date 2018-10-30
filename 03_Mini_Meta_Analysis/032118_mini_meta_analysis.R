# 5-study mini meta-analysis
# Script author: David Chen
# Date: 03/21/18; 09/25/18
# Notes:
# -- Assumption: Sampling variance is small, i.e. 0
# -- Data set had typo, which was fixed on 10/29/18

rm(list=ls());
library(ggplot2);
library(metafor);
library(pheatmap);
MY_COLORS <- colorRampPalette(c("gray88", "gray33", "black"))(128); 
myForestTheme <- theme_bw() +
  theme(axis.text=element_text(size=21,color="black"), 
        axis.title.y=element_blank(), axis.title.x=element_text(size=27,color="black"),
        legend.key.width=unit(2.5,"line"), #legend spacing
        legend.position="top", legend.title=element_text(size=21,color="black",face="bold"), legend.text=element_text(size=18, color="black"));


load5studies <- function(path="~/Dropbox (Christensen Lab)/Christensen Lab - 2018/StatsConsulting_2018/Oncology_Nursing_Meta_Analysis/103018_Studies_included_error_fixed.csv") {
  meta5 <- read.csv(path, row.names=NULL, header=FALSE, stringsAsFactors=FALSE);
  meta5 <- t(meta5);
  colnames(meta5) <- unlist(meta5[1, ]);
  meta5 <- as.data.frame(meta5[-1, ], stringsAsFactors=FALSE);
  meta5$Completed_intervention <- as.numeric(meta5$Completed_intervention);
  meta5$Ratio_of_completed <- as.numeric(meta5$Ratio_of_completed);
  meta5$PA_percent_increase <- as.numeric(meta5$PA_percent_increase);
  meta5$hasGoalSetting <- as.logical(meta5$hasGoalSetting);
  return(meta5);
}

convertToMatForHeat <- function() {
  meta5ForHeat <- load5studies();
  rownames(meta5ForHeat) <- meta5ForHeat$Name;
  meta5ForHeat$Name <- meta5ForHeat$Enrolled <- meta5ForHeat$Completed_intervention <- meta5ForHeat$Region <- NULL;
  meta5ForHeat$Ratio_of_completed <- 100 * meta5ForHeat$Ratio_of_completed; 
  meta5ForHeat$Measurement <- meta5ForHeat$Measurement == "Self-reported";
  for(j in 1:8) {
    meta5ForHeat[ , j] <- 100 * as.logical(meta5ForHeat[ , j]);  
  }
  meta5ForHeat <- as.matrix(t(meta5ForHeat));
  rownames(meta5ForHeat) <- c(
    "Setting: Urban (Yes/No)",
    "Setting: Rural (Yes/No)",
    "Setting: Randomized (Yes/No)",
    "Self-reported (Yes/No)",
    "Final Goal >150 min PA (Yes/No)",
    "Has Goal Setting (Yes/No)",
    "Has BMI Eligibility (Yes/No)",
    "Has PA Eligibility (Yes/No)",
    "Completion Rate (%)",
    "PA Increase (%)"
  );
  return(meta5ForHeat); 
}

runCustomMetaAnalysis <- function(meta5, weightCol="Ratio_of_completed", conditionCol="hasGoalSetting") {
  #'@description Run meta analysis assume variance/heterogeneity is small, i.e. 0
  
  ## All 5 studies:
  print("***********All 5 studies***********")
  res.all5 <- rma(
    yi = meta5$PA_percent_increase,
    vi = rep(0, nrow(meta5)),
    weights = meta5[ , weightCol],
    method = "ML"
  );
  print(res.all5);
  
  ## Subgroup: goal setting
  condition <- meta5[,conditionCol];
  
  print("***********Has goal setting***********");
  res.hasGoalSetting <- rma(
    yi = meta5$PA_percent_increase[condition],
    vi = rep(0, sum(condition)),
    weights = meta5[condition, weightCol],
    method = "ML"
  );
  print(res.hasGoalSetting);
  
  print("***********No goal setting***********");
  res.noGoalSetting <- rma(
    yi = meta5$PA_percent_increase[! condition],
    vi = rep(0, sum(! condition)),
    weights = meta5[! condition, weightCol],
    method = "ML"
  );
  print(res.noGoalSetting);
  
  
  res <- data.frame(
    Category = c("All", "Has Goal Setting", "No Goal Setting"),
    Estimate = c(res.all5$b, res.hasGoalSetting$b, res.noGoalSetting$b),
    ci.lb = c(res.all5$ci.lb, res.hasGoalSetting$ci.lb, res.noGoalSetting$ci.lb),
    ci.ub = c(res.all5$ci.ub, res.hasGoalSetting$ci.ub, res.noGoalSetting$ci.ub),
    P = c(res.all5$pval, res.hasGoalSetting$pval, res.noGoalSetting$pval)
  );
  res$isSig <- res$P < 0.05; 
  
  return(res);
}

main <- function() {
  ## Load data:
  meta5 <- load5studies();
  meta5ForHeat <- convertToMatForHeat();
  
  ## Visualize profiles for all studies:
  png("~/Downloads/Overview_of_5_Studies.png", width=6.5, height=8.27, units="in", res=300);
  pheatmap(
    meta5ForHeat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = MY_COLORS,
    fontsize = 12,
    main = "Overview of Studies"
  );
  dev.off();
  
  ## Run custom meta analysis:
  df4Forest <- runCustomMetaAnalysis(meta5);
  df4Forest$Category <- factor(df4Forest$Category, levels=c("No Goal Setting","Has Goal Setting","All"));
  ggplot(df4Forest, aes(x=Category, y=Estimate, ymin=ci.lb, ymax=ci.ub)) +
    geom_pointrange(size=1.5, position=position_dodge(width=0.5)) +
    geom_hline(yintercept=0, size=0.75, linetype="dashed") +
    labs(y="Avg. % PA increase") +
    scale_color_manual(values=c("royalblue","red")) +
    scale_y_continuous(breaks=seq(-20,100,by=20), limits=c(-20,100)) +
    myForestTheme +
    coord_flip();
  ggsave("~/Downloads/Forest_plot.png", height=8.27, width=11.69);
}

main(); 
