###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TCHA HNSCC covariate data cleaning & model building for multivariable survival analyses including TCGA exp. data 
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")

#------------------------------------------------------------------------------------------
# read data
#------------------------------------------------------------------------------------------

# covariate data 
covs <- readRDS("05_TCGA_expression_survival_analyses/Data_files/TCGA_HNSC_covariates_miRNA_expression_subset.rds")
covs$bcr_patient_uuid <- tolower(covs$bcr_patient_uuid)

#------------------------------------------------------------------------------------------
# clean covariates of interest & explore other variables 
#------------------------------------------------------------------------------------------

# index for laryngeal tumors only 
table(covs$anatomic_neoplasm_subdivision, covs$icd_10)
table(covs$icd_10)

#covs_larynx <- covs[covs$anatomic_neoplasm_subdivision == "Larynx",]
covs_larynx <- covs[covs$icd_10 == "C32.1" | covs$icd_10 == "C32.9",]

# determine if there is enough data to adjust for smoking, drinking and hpv status in data set 
table(is.na(covs_larynx$amount_of_alcohol_consumption_per_day)) 
table(is.na(covs_larynx$number_pack_years_smoked)) 
table(covs_larynx$hpv_status_by_p16_testing) 
# conclusion: too much missingness to adjust for any here 

# consider other variables for inclusion 
table(covs_larynx$postoperative_rx_tx)
table(covs_larynx$hpv_status_by_ish_testing)
table(covs_larynx$drugs)
table(covs_larynx$stage_event)
table(covs_larynx$neoplasm_histologic_grade)
table(covs_larynx$stage_event)
table(covs_larynx$radiation_therapy)

# CONCLUSION: variables for consideration in multivariable cox model: age, sex, race, tumor stage 

#------------------------------------------------------------------------------------------
# model building 
#------------------------------------------------------------------------------------------
# check univariate associatins of variables 
summary(coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status=="Dead") ~ covs_larynx$age_at_initial_pathologic_diagnosis))
# age not strongly associated (will still want to include this, as we know age is often a predictor of survival)
summary(coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status=="Dead") ~ as.factor(covs_larynx$gender)))
# sex is strongly associated 
summary(coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status=="Dead") ~ as.factor(covs_larynx$race_list)))
# race not strongly associated 
summary(coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status=="Dead") ~ as.factor(covs_larynx$tumor_stage)))
# stage not strongly associated 

# make indicator variable numeric so the coxph object made can be used in anova function 
covs_larynx$vital_status_2 <- as.numeric(covs_larynx$vital_status=="Dead")

#~~~~~~~~~~~~~~~~~~
# leave one out model building approach
#~~~~~~~~~~~~~~~~~~

# fit multivariable model 
fit_1 <- coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status_2) ~ 
                 covs_larynx$age_at_initial_pathologic_diagnosis + 
                 covs_larynx$gender + 
                 covs_larynx$race_list + 
                 covs_larynx$tumor_stage)
summary(fit_1)
# sex and race associated 

# remove sex
fit_2 <- coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status_2) ~ 
                 covs_larynx$age_at_initial_pathologic_diagnosis + 
                 #covs_larynx$gender + 
                 covs_larynx$race_list + 
                 covs_larynx$tumor_stage)
summary(fit_2)
# race associated 

# remove race
fit_3 <- coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status_2) ~ 
                 covs_larynx$age_at_initial_pathologic_diagnosis + 
                 as.factor(covs_larynx$gender) + 
                 #as.factor(covs_larynx$race_list) + 
                 as.factor(covs_larynx$tumor_stage))
summary(fit_3)
# sex associated 

fit_4 <- coxph(Surv(covs_larynx$time_to_event, covs_larynx$vital_status_2) ~ 
                 covs_larynx$age_at_initial_pathologic_diagnosis + 
                 as.factor(covs_larynx$gender) + 
                 as.factor(covs_larynx$race_list)) 
#as.factor(covs_larynx$tumor_stage))
summary(fit_4)
# sex associated 


#~~~~~~~~~~~~~~~~~~
# -2loglik testing  
#~~~~~~~~~~~~~~~~~~
# the smaller the value of -2loglik, the better the model explains the observed data 
# if a variable increases -2loglik when ommitted, the variable likely contributes to the outcome
# so use this approach to identify variables that we should keep in the model  

-2*fit_1$loglik[2]
-2*fit_2$loglik[2]
# omition of sex increases -2loglik so reduces model fit, therefore we should keep sex

-2*fit_3$loglik[2]
# omition of race increases -2loglik so reduces model fit, therefore we want to keep race 

-2*fit_4$loglik[2]
# omition of stage increases -2loglik so reduces model fit, therefore we want to keep stage 

# test if there is a significant change in this loglik between models 
anova(fit_2, fit_1) # P = 0.02
anova(fit_3, fit_1) # P = 0.11
anova(fit_4, fit_1) # P = 0.41
# exclusion of sex, race and stage all increase the -2loglik but only the increases in age and race are statistically 
# significant according to liklihood ratio test 

#~~~~~~~~~~~~~~~~~~
# test for correlation between variables  
#~~~~~~~~~~~~~~~~~~
# consider any association between variables

# age vs sex
kruskal.test(covs_larynx$age_at_initial_pathologic_diagnosis ~ 
               as.numeric(covs_larynx$gender == "FEMALE")) # significant 
# age vs race
kruskal.test(covs_larynx$age_at_initial_pathologic_diagnosis ~ 
               as.numeric(covs_larynx$race_list == "NON-WHITE")) # P = 0.42
# age vs stage
kruskal.test(covs_larynx$age_at_initial_pathologic_diagnosis ~ 
               as.numeric(covs_larynx$tumor_stage == "III+IV")) # P = 0.38
# sex vs race
chisq.test(matrix(table(as.numeric(covs_larynx$gender == "FEMALE"), 
                        as.numeric(covs_larynx$race_list == "NON-WHITE")), ncol = 2), simulate.p.value = TRUE)$p.value # non-significant 
# sex vs stage
chisq.test(matrix(table(as.numeric(covs_larynx$gender == "FEMALE"), 
                        as.numeric(covs_larynx$tumor_stage == "III+IV")), ncol = 2), simulate.p.value = TRUE)$p.value # non-significant 
# stage vs race 
chisq.test(matrix(table(as.numeric(covs_larynx$tumor_stage == "III+IV"), 
                        as.numeric(covs_larynx$race_list == "NON-WHITE")), ncol = 2), simulate.p.value = TRUE)$p.value # non-significant 
# note need to use simulate.p.value as values in some cells for chisq tests are <=5 so chisq approx. may be 
# inaccurate so need to apply fisher test

