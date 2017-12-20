
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
# load subvject metadata 
metadata <- read.table("05_TCGA_expression_survival_analyses/Data_files/metadata.cart.2017-04-13T19-36-30.619378.txt", sep = "\t", header = T, stringsAsFactors = F)
# miRNA expression set 
expression_miR <- readRDS("05_TCGA_expression_survival_analyses/Data_files/TCGA_HNSCC_tumor_miRNA_expression_matrix.rds")
# RNA expression set 
expression_mRNA <- readRDS("05_TCGA_expression_survival_analyses/Data_files/TCGA_HNSCC_tumor_RNAseq_expression_matrix.rds")

#------------------------------------------------------------------------------------------
# clean covariates of interest & explore other variables 
#------------------------------------------------------------------------------------------

# confirm subject order still matches 
all(covs$bcr_patient_uuid == rownames(expression_miR))

# add expression of genes from mir-100 locus 
covs$mir_100_continuous <- log2(expression_miR[,which(colnames(expression_miR) == "hsa-mir-100")])
covs$let7a_2_continuous <- log2(expression_miR[,which(colnames(expression_miR) == "hsa-let-7a-2")])
covs$mir125b_1_continuous <- log2(expression_miR[,which(colnames(expression_miR) == "hsa-mir-125b-1")])
# Note: BLID expression must be added separately as expression only available on a subet of subjects 

# index for oral tumors only 
table(covs$anatomic_neoplasm_subdivision)
covs_oral <- covs[covs$anatomic_neoplasm_subdivision == "Oral Cavity",]

# combine death and follow uop variable to get time to event for alive and dead subjects 
covs_oral$time_to_event <- NA
covs_oral$time_to_event[is.na(covs_oral$days_to_death)] <- covs_oral$days_to_last_followup[is.na(covs_oral$days_to_death)]
covs_oral$time_to_event[!is.na(covs_oral$days_to_death)] <- covs_oral$days_to_death[!is.na(covs_oral$days_to_death)]
covs_oral$time_to_event

# consider age, sex, race and tumor stage variables 
covs_oral$age_at_initial_pathologic_diagnosis
covs_oral$gender
covs_oral$gender <- factor(covs_oral$gender, levels = c("MALE", "FEMALE"))
covs_oral$race_list

# process race variable to 'white' 'non-white' 
covs_oral$race_list[covs_oral$race_list == "BLACK OR AFRICAN AMERICAN"] <- "NON-WHITE"
covs_oral$race_list[covs_oral$race_list == "ASIAN"] <- "NON-WHITE"
covs_oral$race_list[covs_oral$race_list == "BLACK OR AFRICAN AMERICAN"] <- "NON-WHITE"
covs_oral$race_list[covs_oral$race_list == ""] <- NA
covs_oral$race_list

# remove those with NAs for race + re-level 
covs_oral$race_list <- factor(covs_oral$race_list, levels = c("WHITE", "NON-WHITE"))
covs_oral_2 <- covs_oral[!is.na(covs_oral$race_list),]

# split complicated stage variable up and extract stage number data 
x <- strsplit(covs_oral_2$stage_event, " ")
stages <- c()
for(i in 1:length(x)){
  stages[i] <- x[[i]][2]
}
stages
table(stages)

# write function to assign each complicated stage character string to stage I, II, III or IV 
tumor_stage <- c()
for(i in 1:length(stages)){
  # if element 2 of character string is I, then if element 3 is I this is stage 3, if else is stage II
  if(strsplit(stages, "")[[i]][2] == "I")
    if(strsplit(stages, "")[[i]][3] == "I"){
      tumor_stage[i] <- "IIIStage"
    } else {tumor_stage[i] <- "IIStage"}
  # if element 2 is V, then stage is 4 
  else if(strsplit(stages, "")[[i]][2] == "V"){
    tumor_stage[i] <- "IVStage"
  } else {tumor_stage[i] <- "IStage"}
}
tumor_stage
table(tumor_stage)

# parse to more clear levels 
tumor_stage[tumor_stage == "IStage"] <- "I+II"
tumor_stage[tumor_stage == "IIStage"] <- "I+II"
tumor_stage[tumor_stage == "IIIStage"] <- "III+IV"
tumor_stage[tumor_stage == "IVStage"] <- "III+IV"

# add tumor stage to covs
covs_oral_2$tumor_stage <- tumor_stage
covs_oral_2$tumor_stage <- factor(covs_oral_2$tumor_stage, levels = c("I+II", "III+IV"))

# check the variables again 
covs_oral_2$gender
covs_oral_2$race_list
covs_oral_2$tumor_stage

# determine if there is enough data to adjust for smoking, drinking and hpv status in larynx data set 
table(is.na(covs_oral_2$amount_of_alcohol_consumption_per_day)) # missing data for 42/72
table(is.na(covs_oral_2$number_pack_years_smoked)) # missing data for 41/72
table(covs_oral_2$hpv_status_by_p16_testing) # missing data for 8/72
# conclusion: too much missingness to adjust for any here 

# consider other variables for inclusion 
table(covs_oral$postoperative_rx_tx)
table(covs_oral$hpv_status_by_ish_testing)
table(covs_oral$drugs)
table(covs_oral$stage_event)
table(covs_oral$neoplasm_histologic_grade)
table(covs_oral$stage_event)
table(covs_oral$radiation_therapy)

# CONCLUSION: variables for consideration in multivariable cox model: age, sex, race, tumor stage 

#------------------------------------------------------------------------------------------
# model building 
#------------------------------------------------------------------------------------------
# check univariate associatins of variables 

summary(coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status=="Dead") ~ covs_oral_2$age_at_initial_pathologic_diagnosis))
# age not strongly associated (will still want to include this, as we know age is often a predictor of survival)
summary(coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status=="Dead") ~ as.factor(covs_oral_2$gender)))
# sex not strongly associated 
summary(coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status=="Dead") ~ as.factor(covs_oral_2$race_list)))
# race is associated 
summary(coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status=="Dead") ~ as.factor(covs_oral_2$tumor_stage)))
# stage not strongly associated (note, even though not associated we may still want to include this as we know stage is associated with survival)

# make indicator variable numeric so the coxph object made can be used in anova function 
covs_oral_2$vital_status_2 <- as.numeric(covs_oral_2$vital_status=="Dead")

#~~~~~~~~~~~~~~~~~~
# leave one out model building approach
#~~~~~~~~~~~~~~~~~~

# fit multivariable model 
fit_1 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 covs_oral_2$gender + 
                 covs_oral_2$race_list + 
                 covs_oral_2$tumor_stage)
summary(fit_1)
# mir_100 expression, age, sex and race associated 
# tumor stage somewhat not associated 

# remove sex
fit_2 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 #covs_oral_2$gender + 
                 covs_oral_2$race_list + 
                 covs_oral_2$tumor_stage)
summary(fit_2)
# mir_100 expression P = 0.077
# race associated 
# age somewhat associated 
# tumor stage not associated 

# remove race
fit_3 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 as.factor(covs_oral_2$gender) + 
                 #as.factor(covs_oral_2$race_list) + 
                 as.factor(covs_oral_2$tumor_stage))
summary(fit_3)
# mir_100 expression P 0.04
# sex almost associated @ 0.05
# age somewhat associated 
# tumor stage not associated 

fit_4 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 as.factor(covs_oral_2$gender) + 
                 as.factor(covs_oral_2$race_list)) 
#as.factor(covs_oral_2$tumor_stage))
summary(fit_4)
# mir_100 expression P 0.03
# age associated 
# race associated 
# sex somewhat associated 

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
anova(fit_2, fit_1) # P = 0.06
anova(fit_3, fit_1) # P = 0.02
anova(fit_4, fit_1) # P = 0.4
# exclusion of sex, race and stage all increase the -2loglik but only the increases in age and race are statistically 
# significant according to liklihood ratio test 

#~~~~~~~~~~~~~~~~~~
# test for correlation between variables  
#~~~~~~~~~~~~~~~~~~
# consider any association between variables

# age vs sex
kruskal.test(covs_oral_2$age_at_initial_pathologic_diagnosis ~ 
               as.numeric(covs_oral_2$gender == "FEMALE")) # significant 
# age vs race
kruskal.test(covs_oral_2$age_at_initial_pathologic_diagnosis ~ 
               as.numeric(covs_oral_2$race_list == "NON-WHITE")) # P = 0.07
# age vs stage
kruskal.test(covs_oral_2$age_at_initial_pathologic_diagnosis ~ 
               as.numeric(covs_oral_2$tumor_stage == "III+IV")) # P = 0.08 
# sex vs race
chisq.test(matrix(table(as.numeric(covs_oral_2$gender == "FEMALE"), 
                        as.numeric(covs_oral_2$race_list == "NON-WHITE")), ncol = 2), simulate.p.value = TRUE)$p.value # non-significant 
# sex vs stage
chisq.test(matrix(table(as.numeric(covs_oral_2$gender == "FEMALE"), 
                        as.numeric(covs_oral_2$tumor_stage == "III+IV")), ncol = 2), simulate.p.value = TRUE)$p.value # non-significant 
# stage vs race 
chisq.test(matrix(table(as.numeric(covs_oral_2$tumor_stage == "III+IV"), 
                        as.numeric(covs_oral_2$race_list == "NON-WHITE")), ncol = 2), simulate.p.value = TRUE)$p.value # non-significant 
# note need to use simulate.p.value as values in some cells for chisq tests are <=5 so chisq approx. may be 
# inaccurate so need to apply fisher test

# CONCLUSION: age, sex, race, & tumor stage should be included in the models 

# save clean covariate data w/ expression included 
saveRDS(covs_oral_2, file = "TCGA_expression_survival_analyses/Data_files/TCGA_oral_cancer_covs_cleaned_w_MIR100HG_expression.rds")
