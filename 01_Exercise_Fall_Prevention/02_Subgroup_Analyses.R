###############################################################################
# RACE data analytics: Table One, Data Visualization, and Subgroup Regressions
# Script author: David (Youdinghuan) Chen
# Dates: 05/18/2017 - 05/21/2017
# Revision: 05/25/2017
###############################################################################

library(coefplot)
library(dplyr)
library(fmsb)
library(gdata)
library(ggplot2)
library(grid)
library(gridExtra)
library(impute)
library(matrixStats)
library(reshape)

setwd("~/Dropbox (Personal)/Dartmouth/ActiveStepRCT_data/")

#------------------------------------Dataset------------------------------------
## Load dataset:
RACE_data      <- read.csv("Race Data.csv")
RACE_data_dict <- read.xls("RACET Data Dictionary.xlsx")

## Check variables in dictionary:
table(RACE_data_dict$Variable.Group.Name)
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Follow-up:  Fall Details"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Demographics"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Follow-up:  Falls"]
RACE_data_dict$Code[RACE_data_dict$Variable.Group.Name=="Baseline Assessment"]

## Motivation:
table(RACE_data$IntervalName)
RACE_data.3mo <- RACE_data[RACE_data$IntervalName=="3-Month Follow-up", ]

sum(!is.na(RACE_data.3mo$times_fallen_3_months))
sum(!is.na(RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode=="ACTIVESTEP"]))
sum(!is.na(RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode=="STANDARD"]))

RACE_data.3mo <- RACE_data.3mo[! is.na(RACE_data.3mo$times_fallen_3_months), , drop=T]
sum(RACE_data.3mo$TreatmentCode=="ACTIVESTEP");
sum(RACE_data.3mo$TreatmentCode=="STANDARD")

## Crude tests:
crude_lm <- lm(times_fallen_3_months ~ factor(TreatmentCode), data=RACE_data.3mo)
summary(crude_lm)
confint(crude_lm)
wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode=="ACTIVESTEP"],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode=="STANDARD"],
  alternative = "less",
  correct = TRUE,
  conf.int = TRUE
)

## This data is only in the screening stage. So match and fill in:
RACE_data.screen <- RACE_data[RACE_data$IntervalName=="Screening", ]
RACE_data.screen <- RACE_data.screen[RACE_data.screen$ReferenceID %in% RACE_data.3mo$ReferenceID, ]

## Important checkpoint: TRUE required
identical(RACE_data.screen$ReferenceID, RACE_data.3mo$ReferenceID) #check: If FALSE, `match` first
identical(match(RACE_data.screen$ReferenceID, RACE_data.3mo$ReferenceID), 1:nrow(RACE_data.3mo)) #double-check

## Update data or covariates missing at 3-month, using data collected at screening:
RACE_data.3mo$peripheral_neuropathy <- RACE_data.screen$peripheral_neuropathy
RACE_data.3mo$stroke                <- RACE_data.screen$stroke
RACE_data.3mo$parkinsons            <- RACE_data.screen$parkinsons
RACE_data.3mo$abc_average_initial   <- RACE_data.screen$abc_average #set initial to ABC at screening
RACE_data.3mo$tug_score_initial     <- RACE_data.screen$tug_score
RACE_data.3mo$dgi_score_initial     <- RACE_data.screen$dgi_score
RACE_data.3mo$visual_impairment     <- RACE_data.screen$visual_impairment
RACE_data.3mo$knee_arthritis        <- RACE_data.screen$knee_arthritis
RACE_data.3mo$hip_arthritis         <- RACE_data.screen$hip_arthritis
RACE_data.3mo$sp_tha                <- RACE_data.screen$sp_tha
RACE_data.3mo$sp_tka                <- RACE_data.screen$sp_tka
RACE_data.3mo$cardiac_arrhythmia    <- RACE_data.screen$cardiac_arrhythmia
RACE_data.3mo$spinal_stenosis       <- RACE_data.screen$spinal_stenosis

## Impute the missing value by kNN approach:
imp_val <- impute.knn(matrix(RACE_data.3mo$abc_average_initial), rowmax=1, colmax=1)$data[is.na(RACE_data.3mo$abc_average_initial)]
RACE_data.3mo$abc_average_initial[is.na(RACE_data.3mo$abc_average_initial)] <- imp_val #replace NA with imputed value
quantile(RACE_data.3mo$abc_average_initial, c(0.6, 0.4))
RACE_data.3mo$abc_average_group <- ifelse(RACE_data.3mo$abc_average_initial >= 66.30, 1, 0)

## Convert to factor level:
RACE_data.3mo$TreatmentCode <- gsub("ACTIVESTEP", 1, RACE_data.3mo$TreatmentCode)
RACE_data.3mo$TreatmentCode <- gsub("STANDARD", 0, RACE_data.3mo$TreatmentCode)
RACE_data.3mo$TreatmentCode <- as.factor(RACE_data.3mo$TreatmentCode)
RACE_data.3mo$parkinsons <- as.factor(RACE_data.3mo$parkinsons)
RACE_data.3mo$peripheral_neuropathy <- as.factor(RACE_data.3mo$peripheral_neuropathy) 
RACE_data.3mo$stroke <- as.factor(RACE_data.3mo$stroke)
RACE_data.3mo$age_group <- as.factor(RACE_data.3mo$age_group)
RACE_data.3mo$abc_average_group <- as.factor(RACE_data.3mo$abc_average_group)

#--------------------------------------------Table One--------------------------------------------
library(tableone)
## Variables to summarize, in the desired order:
vars <- c(
  "sex", "hispa","race", "peripheral_neuropathy","stroke", "parkinsons",
  "visual_impairment", "hip_arthritis", "knee_arthritis", "cardiac_arrhythmia", "spinal_stenosis",
  "age_group", "age_at_screening", "abc_average_initial"
)

## Specify numerically coded variables as factors to stratify
## Not necessary if already factor
factorVars <- c(
  "TreatmentCode", #stratifying variable
  "sex", "hispa", "race",   "age_group",
  "peripheral_neuropathy","stroke", "parkinsons",
  "visual_impairment", "hip_arthritis", "knee_arthritis", "cardiac_arrhythmia", "spinal_stenosis"
) 

## Create tableone object:
tableOne <- CreateTableOne(
  factorVars = factorVars, 
  vars = vars, 
  strata = c("TreatmentCode"), 
  test = FALSE, #do not test for randomized trials
  data = RACE_data.3mo,
  includeNA = TRUE, 
  smd = FALSE
)
tableOne

## Use the quote argument to add quotes to avoid mishandling of spaces. 
## After pasting in Excel, use the TextImportWizard, and press Finish.
# x <- print(tableOne, showAllLevels=TRUE)
# write.csv(x, "~/repos/ActiveStepRCT/LimitedData/052017_ActiveStep3mo.csv", quote=FALSE)

#--------------------------------------------Data visualization--------------------------------------------
## Primary response of interest:
p1 <- ggplot(RACE_data.3mo , aes(TreatmentCode, times_fallen_3_months, color=TreatmentCode)) +
  geom_boxplot(width=0.5) +  #`outlier.shape=NA` to remove outliers
  scale_color_brewer(palette="Set1") + 
  labs(y="Times fallen", title="Num. of times fallen at 3-mo. endpt.") +
  theme_bw()

## Age:
p2 <- ggplot(RACE_data.3mo, aes(TreatmentCode, age_at_screening, color=TreatmentCode)) + 
  geom_boxplot(width=0.5) + 
  scale_color_brewer(palette="Set1") +
  labs(y="Age", title="Age distribution at screening") +
  theme_bw()

## Peripheral neuralpathy proportions: Requires count by dplyr
pn.3mo <- RACE_data.3mo %>%
  select(ReferenceID, TreatmentCode, peripheral_neuropathy) %>%
  group_by(peripheral_neuropathy, TreatmentCode) %>% 
  summarise(n=n())
pn.3mo$peripheral_neuropathy <- gsub(1, "yes", pn.3mo$peripheral_neuropathy)
pn.3mo$peripheral_neuropathy <- gsub(0, "no", pn.3mo$peripheral_neuropathy)
p3 <- ggplot(pn.3mo, aes(x=TreatmentCode, y=n, fill=peripheral_neuropathy)) +
  geom_bar(stat="identity",width=0.5) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_brewer(palette="Accent") +
  labs(x="Presence (1=yes, 0=no)", y="Number of participants", title="Peripheral neuropathy at screening") +
  theme_bw()

## Proportion with PDs by group:
pd.screen <- RACE_data.screen %>%
  select(ReferenceID, TreatmentCode, parkinsons) %>%
  group_by(TreatmentCode, parkinsons) %>% 
  summarise(n=n())
pd.screen$parkinsons <- gsub(1, "yes", pd.screen$parkinsons)
pd.screen$parkinsons <- gsub(0, "no", pd.screen$parkinsons)
p4 <- ggplot(pd.screen, aes(x=TreatmentCode, y=n, fill=parkinsons)) +
  geom_bar(stat="identity",width=0.5) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_brewer(palette="Accent") +
  labs(x="Presence", y="Number of participants", title="Parkinson's Disease at screening") +
  theme_bw()

## Stroke at the time of screening:
stroke.screen <- RACE_data.screen %>%
  select(ReferenceID, TreatmentCode, stroke) %>%
  group_by(TreatmentCode, stroke) %>% 
  summarise(n=n())
stroke.screen$stroke <- gsub(1, "yes", stroke.screen$stroke)
stroke.screen$stroke <- gsub(0, "no", stroke.screen$stroke)
p5 <- ggplot(stroke.screen, aes(x=TreatmentCode, y=n, fill=stroke)) +
  geom_bar(stat="identity",width=0.5) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_brewer(palette="Accent") +
  labs(x="Presence", y="Number of participants", title="Strokes at screening") +
  theme_bw()

pblank <- grid.rect(gp=gpar(col="white")) #placeholder blank plot

grid.arrange(p1,p2,pblank,p3,p4,p5,ncol=3)

#------------------------------------Variable 1: Age Group (>=75 vs. <75)------------------------------------
## Data.frame to store subgroup analysis results:
## Note fmsb::oddsratio() has a slightly different set up of a,b,c,d
res.wilcox <- data.frame(
  Variable  = rep(c("Age 75 or above", "Parkinson's", "Peripheral neu.", "Stroke", "Baseline ABC"), 2),
  Subgroup  = rep(c(0,1), 5),
  test.stat = rep(NA, 10),
  P.value   = rep(NA, 10)
)
res.wilcox <- res.wilcox[order(res.wilcox$Variable), ]
res.em <- NULL

sum(RACE_data.3mo$age_at_screening >= 75); 
sum(RACE_data.3mo$age_group == 1)

## Test
w0_age <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$age_group==0],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0   & RACE_data.3mo$age_group==0],
  alternative = "less",
  correct = TRUE
)
w1_age <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$age_group==1],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0   & RACE_data.3mo$age_group==1],
  alternative = "less",
  correct = TRUE
)

res.wilcox$test.stat[res.wilcox$Variable=="Age 75 or above" & res.wilcox$Subgroup==0] <- w0_age$statistic
res.wilcox$P.value[res.wilcox$Variable=="Age 75 or above" & res.wilcox$Subgroup==0]   <- w0_age$p.value
res.wilcox$test.stat[res.wilcox$Variable=="Age 75 or above" & res.wilcox$Subgroup==1] <- w1_age$statistic
res.wilcox$P.value[res.wilcox$Variable=="Age 75 or above" & res.wilcox$Subgroup==1]   <- w1_age$p.value

em.age <- lm(times_fallen_3_months ~ age_group*TreatmentCode, data=RACE_data.3mo)
summary(em.age)$coef
confint(em.age)

#------------------------------------Variable 2: Peripheral Neuropathy------------------------------------
w1_pn <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$peripheral_neuropathy==1],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0 & RACE_data.3mo$peripheral_neuropathy==1],
  alternative = "less",
  correct = TRUE
)
w0_pn <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$peripheral_neuropathy==0],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0 & RACE_data.3mo$peripheral_neuropathy==0],
  alternative = "less",
  correct = TRUE
)

## Updates:
res.wilcox$test.stat[res.wilcox$Variable=="Peripheral neu." & res.wilcox$Subgroup==0] <- w0_pn$statistic
res.wilcox$P.value[res.wilcox$Variable=="Peripheral neu." & res.wilcox$Subgroup==0]   <- w0_pn$p.value
res.wilcox$test.stat[res.wilcox$Variable=="Peripheral neu." & res.wilcox$Subgroup==1] <- w1_pn$statistic
res.wilcox$P.value[res.wilcox$Variable=="Peripheral neu." & res.wilcox$Subgroup==1]   <- w1_pn$p.value

em.pn <- lm(times_fallen_3_months ~ peripheral_neuropathy*TreatmentCode, data=RACE_data.3mo)
summary(em.pn)
confint(em.pn)

#------------------------------------Variable 3: stroke------------------------------------
w1_stroke <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$stroke==1],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0 & RACE_data.3mo$stroke==1],
  alternative = "less",
  correct = TRUE
)
w0_stroke <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$stroke==0],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0 & RACE_data.3mo$stroke==0],
  alternative = "less",
  correct = TRUE
)

## Updates:
res.wilcox$test.stat[res.wilcox$Variable=="Stroke" & res.wilcox$Subgroup==0] <- w0_stroke$statistic
res.wilcox$P.value[res.wilcox$Variable=="Stroke" & res.wilcox$Subgroup==0]   <- w0_stroke$p.value
res.wilcox$test.stat[res.wilcox$Variable=="Stroke" & res.wilcox$Subgroup==1] <- w1_stroke$statistic
res.wilcox$P.value[res.wilcox$Variable=="Stroke" & res.wilcox$Subgroup==1]   <- w1_stroke$p.value

em.stroke <- lm(times_fallen_3_months ~ stroke*TreatmentCode, data=RACE_data.3mo)
summary(em.stroke)$coef
confint(em.stroke)

#------------------------------------Variable 4: Parkinson's disease------------------------------------
w1_parkin <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$parkinsons==1],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0   & RACE_data.3mo$parkinsons==1],
  alternative = "less",
  correct = TRUE
)
w0_parkin <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$parkinsons==0],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0 & RACE_data.3mo$parkinsons==0],
  alternative = "less",
  correct = TRUE
)

## Updates:
res.wilcox$test.stat[res.wilcox$Variable=="Parkinson's" & res.wilcox$Subgroup==0] <- w0_parkin$statistic
res.wilcox$P.value[res.wilcox$Variable=="Parkinson's" & res.wilcox$Subgroup==0]   <- w0_parkin$p.value
res.wilcox$test.stat[res.wilcox$Variable=="Parkinson's" & res.wilcox$Subgroup==1] <- w1_parkin$statistic
res.wilcox$P.value[res.wilcox$Variable=="Parkinson's" & res.wilcox$Subgroup==1]   <- w1_parkin$p.value

em.parkin <- lm(times_fallen_3_months ~ parkinsons*TreatmentCode, data=RACE_data.3mo)
summary(em.parkin)$coef
confint(em.parkin)

#------------------------------------Variable 5: Baseline ABC group------------------------------------
w1_abc <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$abc_average_group==1],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0   & RACE_data.3mo$abc_average_group==1],
  alternative = "less",
  correct = TRUE
)
w0_abc <- wilcox.test(
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==1 & RACE_data.3mo$abc_average_group==0],
  RACE_data.3mo$times_fallen_3_months[RACE_data.3mo$TreatmentCode==0   & RACE_data.3mo$abc_average_group==0],
  alternative = "less",
  correct = TRUE
)

## Updates:
res.wilcox$test.stat[res.wilcox$Variable=="Baseline ABC" & res.wilcox$Subgroup==0] <- w0_abc$statistic
res.wilcox$P.value[res.wilcox$Variable=="Baseline ABC" & res.wilcox$Subgroup==0] <- w0_abc$p.value
res.wilcox$test.stat[res.wilcox$Variable=="Baseline ABC" & res.wilcox$Subgroup==1] <- w1_abc$statistic
res.wilcox$P.value[res.wilcox$Variable=="Baseline ABC" & res.wilcox$Subgroup==1] <- w1_abc$p.value

em.abc <- lm(times_fallen_3_months ~ abc_average_group*TreatmentCode, data=RACE_data.3mo)
summary(em.abc)$coef
confint(em.abc)

#------------------------------------Data visualization------------------------------------
res.em <- as.data.frame(rbind(
  cbind(summary(em.age)   $coef[3, ,drop=F], confint(em.age)[3, ,drop=F]),
  cbind(summary(em.parkin)$coef[3, ,drop=F], confint(em.parkin)[3, ,drop=F]),
  cbind(summary(em.pn)    $coef[3, ,drop=F], confint(em.pn)[3, ,drop=F]),
  cbind(summary(em.stroke)$coef[3, ,drop=F], confint(em.stroke)[3, ,drop=F]), 
  cbind(summary(em.abc)   $coef[3, ,drop=F], confint(em.abc)[3, ,drop=F])
), row.names=1:5)
rownames(res.em) <- c("age", "Parkinson's", "Peripheral neuropathy", "stroke", "ABC")
colnames(res.em)

em.plt <- ggplot(res.em, aes(x=rownames(res.em), y=Estimate, ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2, color="darkgray") +
  labs(x="Interaction",y="Coef. estimate (95% CI)", title="Effect Modification by Subgroup") + 
  coord_flip() +
  theme_bw()

res.em$`Pr(>|t|)`<- round(res.em$`Pr(>|t|)`, 3)
res.em$Estimate  <- round(res.em$Estimate,   2)
res.em$`2.5 %`   <- round(res.em$`2.5 %`,    2)
res.em$`97.5 %`  <- round(res.em$`97.5 %`,   2)
res.table <- data.frame(
  Interaction   = rownames(res.em),
  `Coefficient` = res.em$Estimate,
  `P-value`     = res.em$`Pr(>|t|)`,
  CI.lower      = res.em$`2.5 %`,
  CI.upper      = res.em$`97.5 %`
)
cex <-  0.75
grobTheme <- ttheme_default(
  core    = list(fg_params=list(cex=cex), bg_params = list(fill=NA, col=1)),
  colhead = list(fg_params=list(cex=cex), bg_params = list(fill=NA, col=1)),
  rowhead = list(fg_params=list(cex=cex), bg_params = list(fill=NA, col=1))
)
em.tbl <- tableGrob(res.table, rows=NULL, theme=grobTheme)
grid.arrange(em.plt, em.tbl, ncol=2)

res.wilcox[ , -1] <- round(res.wilcox[ , -1], 2)
grid.arrange(tableGrob(res.wilcox, theme=grobTheme, rows=NULL))

#------------------------------------Linear & generalized linear model, fixed effect------------------------------------
library(sandwich)
library(lmtest)

## Linear model with adjustment of relevant covariate
lmod5 <- lm(times_fallen_3_months~TreatmentCode+age_group+parkinsons+peripheral_neuropathy+stroke+abc_average_group, data=RACE_data.3mo)
svar.lm <- sandwich(lmod5)
res.lmod5 <- as.data.frame(cbind(
  coeftest(lmod5, vcov.=svar.lm)[2:7, ],
  coefci(lmod5, vcov.=svar.lm)[2:7, ]
))
c0 <- ggplot(res.lmod5, aes(x=rownames(res.PoisMod), y=Estimate, ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange(color="mediumorchid") + 
  geom_hline(yintercept=0, lty=2, color="darkgray") +
  labs(x="Variable",y="Coef. estimate (95% CI)", title='Generalized linear model (GLM) with robust var') + 
  coord_flip() +
  theme_bw()
res.lmod5 <- round(res.lmod5, 2)
c0.tab <- tableGrob(res.lmod5, theme=grobTheme)

poisMod5 <- glm(
  times_fallen_3_months ~ TreatmentCode+age_group+parkinsons+peripheral_neuropathy+stroke+abc_average_group,
  data = RACE_data.3mo,
  family = poisson
)
svar.Pois <- sandwich(poisMod5) #sandwich variance
res.PoisMod <- as.data.frame(cbind(
  coeftest(poisMod5, vcov.=svar.Pois)[2:7, ],
  coefci(poisMod5, vcov.=svar.Pois)[2:7, ]
))
c2 <- ggplot(res.PoisMod, aes(x=rownames(res.PoisMod), y=Estimate, ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange(color="darkolivegreen3") + 
  geom_hline(yintercept=0, lty=2, color="darkgray") +
  labs(x="Variable",y="Coef. estimate (95% CI)", title='Poisson model with robust var') + 
  coord_flip() +
  theme_bw()
res.PoisMod <- round(res.PoisMod, 2)
c2.tab <- tableGrob(res.PoisMod, theme=grobTheme)

grid.arrange(c0, c0.tab, c2, c2.tab, ncol=2)

combined.coefs <- data.frame(
  `linear model coefficient`  = summary(lmod5)$coef[,1],
  `linear model P`            = summary(lmod5)$coef[,4],
  `Poisson model coefficient` = summary(poisMod5)$coef[,1],
  `Poisson model P`           = summary(poisMod5)$coef[,4]
)
combined.coefs <- round(combined.coefs, 3)
# write.csv(combined.coefs, file="~/Dropbox (Personal)/QBS121_Assignments/QBS121_term_proj/coefficients_from_adjusted_models.csv",quote=F,row.names=T)

