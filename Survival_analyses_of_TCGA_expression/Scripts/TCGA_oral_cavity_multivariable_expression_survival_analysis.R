
###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test association between survival time in HNSCC TCGA subjects and tumnor tissue gene expression of all 4 genes 
# (mir100, mir1251b, let7a2, BLID) derived from the long non-coding RNA locus MIR100HG
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")
library(survival)
library(survminer)

#------------------------------------------------------------------------------------------
# load data
#------------------------------------------------------------------------------------------
covs <- readRDS("TCGA_expression_survival_analyses/Data_files/TCGA_oral_cancer_covs_cleaned_w_MIR100HG_expression.rds")

# covariate file for BLID contains a different no of subjects as this is the subset of subjects with available RNA-seq data 
covs_RNA <- readRDS("TCGA_expression_survival_analyses/Data_files/TCGA_oral_cancer_covs_cleaned_w_BLID_expression.rds")

#------------------------------------------------------------------------------------------
# fit multivariable cox models for each gene at MIR100HG locus 
#------------------------------------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~
# mir-100
#~~~~~~~~~~~~~~~~~~

# fit model w/ mir_100 as continuous variable 
fit_cont <- coxph(Surv(covs$time_to_event, covs$vital_status_2) ~ 
                    covs$mir_100_continuous 
                  + covs$age_at_initial_pathologic_diagnosis 
                  + covs$gender 
                  + covs$race_list
                  + covs$tumor_stage)
summary(fit_cont)

# build results table 
results <- matrix(NA, nrow = 4, ncol = 8)
colnames(results) <- c("Un_P", "Un_HR", "Un_95_CI_LB", "Un_95_CI_UB", 
                              "Ad_P", "Ad_HR", "Ad_95_CI_LB", "Ad_95_CI_UB")
rownames(results) <- c("mir100", "let7a2", "mir125b1", "blid")

# fill for adjusted model 
results[1, "Ad_P"] <- round(summary(fit_cont)$coef[1,5], digits = 8)
results[1, "Ad_HR"] <- round(exp(summary(fit_cont)$coef[1,1]), digits = 8)
results[1, "Ad_95_CI_LB"] <- round(exp(confint.default(fit_cont))[1, ][1], digits = 8)
results[1, "Ad_95_CI_UB"] <- round(exp(confint.default(fit_cont))[1, ][2], digits = 8)

# fit unadjusted model w/ mir_100 as categorial variable 
fit_con_unadj <- coxph(Surv(covs$time_to_event, covs$vital_status_2) ~ 
                         covs$mir_100_continuous)

# fill for unadjusted model
results[1, "Un_P"] <- round(summary(fit_con_unadj)$coef[1,5], digits = 8)
results[1, "Un_HR"] <- round(exp(summary(fit_con_unadj)$coef[1,1]), digits = 8)
results[1, "Un_95_CI_LB"] <- round(exp(confint.default(fit_con_unadj))[1, ][1], digits = 8)
results[1, "Un_95_CI_UB"] <- round(exp(confint.default(fit_con_unadj))[1, ][2], digits = 8)
rm(fit_cont)

#~~~~~~~~~~~~~~~~~~
# let7a2
#~~~~~~~~~~~~~~~~~~

# fit model w/ mir_100 as continuous variable 
fit_cont <- coxph(Surv(covs$time_to_event, covs$vital_status_2) ~ 
                    covs$let7a_2_continuous 
                  + covs$age_at_initial_pathologic_diagnosis 
                  + covs$gender 
                  + covs$race_list
                  + covs$tumor_stage)
summary(fit_cont)

# fill for adjusted model 
results[2, "Ad_P"] <- round(summary(fit_cont)$coef[1,5], digits = 8)
results[2, "Ad_HR"] <- round(exp(summary(fit_cont)$coef[1,1]), digits = 8)
results[2, "Ad_95_CI_LB"] <- round(exp(confint.default(fit_cont))[1, ][1], digits = 8)
results[2, "Ad_95_CI_UB"] <- round(exp(confint.default(fit_cont))[1, ][2], digits = 8)

# fit unadjusted model w/ mir_100 as categorial variable 
fit_con_unadj <- coxph(Surv(covs$time_to_event, covs$vital_status_2) ~ 
                         covs$let7a_2_continuous)

# fill for unadjusted model
results[2, "Un_P"] <- round(summary(fit_con_unadj)$coef[1,5], digits = 8)
results[2, "Un_HR"] <- round(exp(summary(fit_con_unadj)$coef[1,1]), digits = 8)
results[2, "Un_95_CI_LB"] <- round(exp(confint.default(fit_con_unadj))[1, ][1], digits = 8)
results[2, "Un_95_CI_UB"] <- round(exp(confint.default(fit_con_unadj))[1, ][2], digits = 8)
rm(fit_cont)

#~~~~~~~~~~~~~~~~~~
# mir-125b-1
#~~~~~~~~~~~~~~~~~~

# fit model w/ mir_100 as continuous variable 
fit_cont <- coxph(Surv(covs$time_to_event, covs$vital_status_2) ~ 
                    covs$mir125b_1_continuous 
                  + covs$age_at_initial_pathologic_diagnosis 
                  + covs$gender 
                  + covs$race_list
                  + covs$tumor_stage)
summary(fit_cont)

# fill for adjusted model 
results[3, "Ad_P"] <- round(summary(fit_cont)$coef[1,5], digits = 8)
results[3, "Ad_HR"] <- round(exp(summary(fit_cont)$coef[1,1]), digits = 8)
results[3, "Ad_95_CI_LB"] <- round(exp(confint.default(fit_cont))[1, ][1], digits = 8)
results[3, "Ad_95_CI_UB"] <- round(exp(confint.default(fit_cont))[1, ][2], digits = 8)

# fit unadjusted model w/ mir_100 as categorial variable 
fit_con_unadj <- coxph(Surv(covs$time_to_event, covs$vital_status_2) ~ 
                         covs$mir125b_1_continuous)

# fill for unadjusted model
results[3, "Un_P"] <- round(summary(fit_con_unadj)$coef[1,5], digits = 8)
results[3, "Un_HR"] <- round(exp(summary(fit_con_unadj)$coef[1,1]), digits = 8)
results[3, "Un_95_CI_LB"] <- round(exp(confint.default(fit_con_unadj))[1, ][1], digits = 8)
results[3, "Un_95_CI_UB"] <- round(exp(confint.default(fit_con_unadj))[1, ][2], digits = 8)
rm(fit_cont)

#~~~~~~~~~~~~~~~~~~
# BLID
#~~~~~~~~~~~~~~~~~~
covs_RNA$vital_status_2 <- as.numeric(covs_RNA$vital_status=="Dead")
# fit model w/ mir_100 as continuous variable 
fit_cont <- coxph(Surv(covs_RNA$time_to_event, covs_RNA$vital_status_2) ~ 
                    covs_RNA$blid_continuous 
                  + covs_RNA$age_at_initial_pathologic_diagnosis 
                  + covs_RNA$gender 
                  + covs_RNA$race_list
                  + covs_RNA$tumor_stage)
summary(fit_cont)

# fill for adjusted model 
results[4, "Ad_P"] <- round(summary(fit_cont)$coef[1,5], digits = 8)
results[4, "Ad_HR"] <- round(exp(summary(fit_cont)$coef[1,1]), digits = 8)
results[4, "Ad_95_CI_LB"] <- round(exp(confint.default(fit_cont))[1, ][1], digits = 8)
results[4, "Ad_95_CI_UB"] <- round(exp(confint.default(fit_cont))[1, ][2], digits = 8)

# fit unadjusted model w/ mir_100 as categorial variable 
fit_con_unadj <- coxph(Surv(covs_RNA$time_to_event, covs_RNA$vital_status_2) ~ 
                         covs_RNA$blid_continuous)
summary(fit_con_unadj)

# fill for unadjusted model
results[4, "Un_P"] <- round(summary(fit_con_unadj)$coef[1,5], digits = 8)
results[4, "Un_HR"] <- round(exp(summary(fit_con_unadj)$coef[1,1]), digits = 8)
results[4, "Un_95_CI_LB"] <- round(exp(confint.default(fit_con_unadj))[1, ][1], digits = 8)
results[4, "Un_95_CI_UB"] <- round(exp(confint.default(fit_con_unadj))[1, ][2], digits = 8)

# save results 
write.csv(results, 'TCGA_expression_survival_analyses/Data_files/TCGA_coxph_model_results_expression.csv')


