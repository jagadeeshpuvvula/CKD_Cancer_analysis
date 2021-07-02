#load libs
library(survival)
library(survminer)
library(dplyr)
library(arules)

#Load dataset
ckd<- read.csv("C:/Users/jagad/Desktop/CKD_CA_final/Final_analysis/CKDCA_DB.csv",
               sep= "," , header = TRUE)
ckd$gfr_cat <- cut(ckd$GFR_F, c(0,29,44,59,134)) #group categories adapted from literature
#table((discretize(ckd$GFR_F, method = "cluster", breaks = 5)))
ckd$age_grp <- cut(ckd$age_at_ckd_dx, c(20,59,69,79, 105))
str(ckd)

summary(ckd)

surv_object <- Surv(time = ckd$futime_ca, event = ckd$fstat_ca)
surv_object

##################################################################################
#time to cancer by GFR group
fit1 <- survfit(surv_object ~ gfr_cat, data = ckd)
summary(fit1)
ggsurvplot(fit1, data = ckd, pval = TRUE, title="Time to develop cancer among CKD patients", 
           legend.title="GFR groups", xscale="d_y", risk.table = T, surv.median.line = "hv")

# Fit a Cox proportional hazards model: time to cancer by GFR
ckd$gfr_cat = relevel(ckd$gfr_cat, ref = "(0,29]")
fit.coxph <- coxph(surv_object ~ gfr_cat, data = ckd)
ggforest(fit.coxph, data = ckd)
##################################################################################

#Time to cancer by gender
fit1 <- survfit(surv_object ~ SEX_CD, data = ckd)
summary(fit1)
ggsurvplot(fit1, data = ckd, pval = TRUE)

# Fit a Cox proportional hazards model: time to cancer by GFR
ckd$SEX_CD = relevel(ckd$SEX_CD, ref = "Female")
fit.coxph <- coxph(surv_object ~ SEX_CD, data = ckd)
ggforest(fit.coxph, data = ckd)
##################################################################################


#Time to cancer by race
fit2 <- survfit(surv_object ~ RACE_CD, data = ckd)
summary(fit1)
ggsurvplot(fit2, data = ckd, pval = TRUE)

# Fit a Cox proportional hazards model
ckd$RACE_CD = relevel(ckd$RACE_CD, ref = "White or Caucasian")
fit.coxph <- coxph(surv_object ~ RACE_CD, data = ckd)
ggforest(fit.coxph, data = ckd)

###################################################################################
###################################################################################

############# LUNG CANCER ################
lungca<- filter(ckd, ICD_Group == "C34" | ICD_Group =="NO_Cancer")
write.csv(lungca, file = "C:/Users/jagad/Desktop/CKD_CA_final/Final_working_files/CKD_cohort_v4.csv"
          , row.names = FALSE)

surv_object <- Surv(time = lungca$�..futime_death, event = lungca$fustat_death)
surv_object

#survival time in CKD patients with and without cancer diagnosis
fit1 <- survfit(surv_object ~ ICD_Group, data = lungca)
summary(fit1)
ggsurvplot(fit1, data = lungca, pval = TRUE)

# Fit a Cox proportional hazards model
lungca$ICD_Group = relevel(lungca$ICD_Group, ref = "NO_Cancer")
fit.coxph <- coxph(surv_object ~ ICD_Group$C34, data = lungca)
ggforest(fit.coxph, data = lungca) 
