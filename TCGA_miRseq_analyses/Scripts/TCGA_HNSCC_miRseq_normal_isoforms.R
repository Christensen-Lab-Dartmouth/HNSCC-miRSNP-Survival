
#######################

# Analyze miRNAseq data from TCGA for only normal adjacent tissue samples 

#######################

# Notes: Data downloaded from GDC website using the 'sample type' filter for: 
# primary tumor +
# miRNAseq 
# TCGA_HNSC 

#######################
# Authors: Owen Wilkins 

rm(list = ls())
setwd("/Users/Owen/Downloads/")
library(survival)
library(survminer)

#------------------------------------------------------------------------------------------
# Load data
#------------------------------------------------------------------------------------------

covs <- read.csv('/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/7_TCGA_expression/Data_files/TCGA_HNSC_covariates.csv', header = T, stringsAsFactors = F)
covs$bcr_patient_uuid <- tolower(covs$bcr_patient_uuid)

setwd("/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/5_TCGA_miRseq/Data_files/TCGA_HNSC_miRNA-seq/")
dir("/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/5_TCGA_miRseq/Data_files/TCGA_HNSC_miRNA-seq/")

metadata <- read.table("/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/5_TCGA_miRseq/Data_files/metadata.cart.2017-04-13T19-36-30.619378.txt", sep = "\t", header = T, stringsAsFactors = F)

# get all file name in directory 
filelist = list.files(path = dir(),pattern = ".txt", full.names = TRUE)
filelist

# note: there are "annotation.txt" files associated with 36 subjects that indicate if case was found to be a recurrence
# after submission, if the patient had a simulatneous cancer of another kind etc. 
# these cases should be removed 
indicies_to_drop <- grep("annotations.txt", filelist)
filelist_2 <- filelist[-indicies_to_drop]

# read in all files 
datalist = lapply(filelist_2, function(x)read.table(x, header=T, stringsAsFactors = F)) 

# split the IDs from filelist_2 and match them to the metadata 
get_ids <- function(vector_to_split) strsplit(vector_to_split, "/")[[1]][1]
filelist_3 <- sapply(filelist_2, get_ids)
names(datalist) <- filelist_3
names(datalist)[1:5]
names(filelist_3) <- NULL

# check that all filelist items are in metadata 
table(filelist_3 %in% metadata$cases.case_id)
# 1 is missing as filelist has 501 eklements and metadata data has 500

# check if file_id variable in metadata corresponds to file identifier for TCGA miRSeq data files 
metadata$file_id %in% filelist_3

# check if the cases.case_id variable in metadata is the patient uuid 
metadata$cases.case_id %in% covs$bcr_patient_uuid
# yes it is 

#filelist_4 <- filelist_3[match(metadata_2$file_name, filelist_3)] # index filelist for subjects w/ metadata 
filelist_4 <- filelist_3[order(filelist_3)] # order metadata and filelist the same way  
metadata_2 <- metadata[order(metadata$file_id),]
all(metadata_2$file_id == filelist_4) # check they are identical 

#datalist_2 <- datalist[match(filelist_4, names(datalist))] # index datalist for elements w/ matching IDs in filelist_5
datalist_2 <- datalist[order(names(datalist))] # order datalist_2 same way as filelist_5
all(filelist_4 == names(datalist_2)) # check they are identical 

names(datalist_2)[1:5]
filelist_4[1:5]
metadata_2$file_id[1:5]

names(datalist_2) <- metadata_2$cases.case_id # give datalist patient_uuid as names (instead of file uuid)
datalist_3 <- datalist_2[order(names(datalist_2))] # order based on this new variable 
names(datalist_3) %in% covs$bcr_patient_uuid

# index and order covs for only those subjects in datalist_3 
covs_2 <- covs[match(names(datalist_3), covs$bcr_patient_uuid),]
covs_3 <- covs_2[order(covs_2$bcr_patient_uuid),]

all(covs_3$bcr_patient_uuid == names(datalist_3))

expression <- matrix(NA, nrow = length(covs_3$bcr_patient_uuid), ncol = length(datalist_3[[1]]$miRNA_ID)) # make matrix to hold expression data 
colnames(expression) <- datalist_3[[1]]$miRNA_ID
rownames(expression) <- covs_2$bcr_patient_uuid
for(i in 1:length(datalist_3)){
  expression[i,] <- datalist_3[[i]]$reads_per_million_miRNA_mapped
} # fill matrix w/ expression data 
expression[1:10, 1:5]

# clean up workspace
rm(filelist_4, filelist_3, filelist_2, filelist, covs_2, covs, metadata, metadata_2, 
   datalist, datalist_2, datalist_3, indicies_to_drop)

#------------------------------------------------------------------------------------------
# extract and save average reads per million (rpm) for each miR of interest 
#------------------------------------------------------------------------------------------

# check expression of miRs of interest 
median(log2(expression[, "hsa-mir-100"]))
median(log2(expression[, "hsa-mir-125b-1"]))
median(log2(expression[, "hsa-let-7a-2"]))

hist(log2(expression[, "hsa-mir-100"]))
hist(log2(expression[, "hsa-mir-125b-1"]))
hist(log2(expression[, "hsa-let-7a-2"]))

plot(log2(expression[, "hsa-mir-100"]), log2(expression[, "hsa-mir-125b-1"]))
plot(log2(expression[, "hsa-mir-100"]), log2(expression[, "hsa-let-7a-2"]))
plot(log2(expression[, "hsa-mir-125b-1"]), log2(expression[, "hsa-let-7a-2"]))

summary(lm(log2(expression[, "hsa-mir-100"]) ~ log2(expression[, "hsa-mir-125b-1"])))
summary(lm(log2(expression[, "hsa-mir-100"]) ~ log2(expression[, "hsa-let-7a-2"])))
summary(lm(log2(expression[, "hsa-mir-125b-1"]) ~ log2(expression[, "hsa-let-7a-2"])))

cor(log2(expression[, "hsa-mir-100"]), log2(expression[, "hsa-mir-125b-1"]))
cor(log2(expression[, "hsa-mir-100"]), log2(expression[, "hsa-let-7a-2"]))
cor(log2(expression[, "hsa-mir-125b-1"]), log2(expression[, "hsa-let-7a-2"]))


#
mir100 <- data.frame(log2(expression[, "hsa-mir-100"]))
mir125b1 <- data.frame(log2(expression[, "hsa-mir-125b-1"]))
mirlet7a2 <- data.frame(log2(expression[, "hsa-let-7a-2"]))

mir100$miRNA <- "miR-100"
colnames(mir100) <- c("expression", "miRNA")
mir125b1$miRNA <- "miR-125b-1"
colnames(mir125b1) <- c("expression", "miRNA")
mirlet7a2$miRNA <- "miR-let-7a-2"
colnames(mirlet7a2) <- c("expression", "miRNA")

mirs_exp <- rbind(mir100, mir125b1, mirlet7a2)


ggplot(mirs_exp, aes(expression, fill = miRNA)) + geom_density(alpha = 0.2)


carrots <- data.frame(length = rnorm(100000, 6, 2))
cukes <- data.frame(length = rnorm(50000, 7, 2.5))

#Now, combine your two dataframes into one.  First make a new column in each that will be a variable to identify where they came from later.
carrots$veg <- 'carrot'
cukes$veg <- 'cuke'

#and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes)
ggplot(vegLengths, aes(length, fill = veg)) + geom_density(alpha = 0.2)


#~~~~~~~~~~~~~~~ mir_100 ~~~~~~~~~~~~~~~#
mir_100 <- log2(expression[,which(colnames(expression) == "hsa-mir-100")])
covs_3$mir_100_continuous <- mir_100 # add continuous expressioin variable to covs 
plot(mir_100)

# EAF of rs3777710 in mir_100 was 0.08 so should divide subjects up into expression groups that reflect this 
# given this MAF of 0.08, could divide subjects into 12 groups based on expression of mir_100 
# calculate percentiles for mir_100 expression 
perc.rank <- ecdf(mir_100)
percentile <- sapply(mir_100, perc.rank)

# index top 8% and bottom 8% of mir_100 expressors 
mir_100_1 <- percentile[percentile >= (1-0.5)]
mir_100_2 <- percentile[percentile >= (1-(1)) & percentile < (1-0.5)]

# check to make sure expression values divided correctly 
boxplot(mir_100_1, mir_100_2)

# add indicator variable to covs to assign subjects to mir_100 expression groups 
covs_3$mir_100 <- NA
covs_3$mir_100[covs_3$bcr_patient_uuid %in% names(mir_100_1)]  <- "G1"
covs_3$mir_100[covs_3$bcr_patient_uuid %in% names(mir_100_2)]  <- "G2"
covs_3$mir_100

# index for laryngeal tumors only 
table(covs_3$anatomic_neoplasm_subdivision)
covs_oral <- covs_3[covs_3$anatomic_neoplasm_subdivision == "Oral Cavity",]

# combine death and follow uop variable to get time to event for alive and dead subjects 
covs_oral$time_to_event <- NA
covs_oral$time_to_event[is.na(covs_oral$days_to_death)] <- covs_oral$days_to_last_followup[is.na(covs_oral$days_to_death)]
covs_oral$time_to_event[!is.na(covs_oral$days_to_death)] <- covs_oral$days_to_death[!is.na(covs_oral$days_to_death)]
covs_oral$time_to_event

# generate Kaplan-Meier object for mir_100 
fit <- survfit(Surv(covs_oral$time_to_event, covs_oral$vital_status=="Dead") ~ as.factor(covs_oral$mir_100))

ggsurvplot(fit, data = covs_oral,
           title = "", 
           legend.title = "",
           xlab = "Time (days)",
           ylab = "Survival probability",
           #xlim = c(0, 260),
           # Change legends: title & labels
           #legend.labs = c("mir_100 High", "mir_100 Low"),
           # Add p-value and confidence intervals
           pval = TRUE,
           conf.int = FALSE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           #palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_bw())


#------------------------------------------------------------------------------------------
# prepare variables for cox model analysis 
#------------------------------------------------------------------------------------------

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

# CONCLUSION: variables for consideration in multivariate cox model: age, sex, race, tumor stage 

#------------------------------------------------------------------------------------------
# model building 
#------------------------------------------------------------------------------------------

# consider univariate associatioins for potential inclusion in model 
summary(coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status=="Dead") ~ covs_oral_2$mir_100))
# mir_100 expression as a categorical variable somewhat associated 
summary(coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status=="Dead") ~ covs_oral_2$mir_100_continuous))
# mir_100 expression as a continuous variable somewhat associated 
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

# fit multivariate model 
fit_1 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$mir_100 + 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 covs_oral_2$gender + 
                 covs_oral_2$race_list + 
                 covs_oral_2$tumor_stage)
summary(fit_1)

# mir_100 expression, age, sex and race associated 
# tumor stage somewhat not associated 

# remove sex
fit_2 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$mir_100 + 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 #covs_oral_2$gender + 
                 covs_oral_2$race_list + 
                 covs_oral_2$tumor_stage)
summary(fit_2)
# mir_100 expression P = 0.16
# race associated 
# age somewhat associated 
# tumor stage not associated 

# remove race
fit_3 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$mir_100 + 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 as.factor(covs_oral_2$gender) + 
                 #as.factor(covs_oral_2$race_list) + 
                 as.factor(covs_oral_2$tumor_stage))
summary(fit_3)
# mir_100 expression P 0.08
# sex almost associated @ 0.05
# age somewhat associated 
# tumor stage not associated 

fit_4 <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                 covs_oral_2$mir_100 + 
                 covs_oral_2$age_at_initial_pathologic_diagnosis + 
                 as.factor(covs_oral_2$gender) + 
                 as.factor(covs_oral_2$race_list)) 
#as.factor(covs_oral_2$tumor_stage))
summary(fit_4)
# mir_100 expression P 0.09
# age somewhat associated 
# race associated 
# sex somewhat associated 

# the smaller the value of -2loglik, the better the model explains the observed data 
# so if something increases the -loglik when omitted, that variable is actually something we probably want top keep in the model 
# so we can discard variables that do not significantly increase -2loglik when omitted, because they are probably not important for model fit
# alternatively, if something decreases -2loglik when omitted, then that variable is probably not important 
# if when you include a varibale in a model, the -2loglik decreases, it is likely something you want to include 

-2*fit_1$loglik[2]
-2*fit_2$loglik[2]
# omition of sex increases -2loglik so reduces model fit, therefore we want to keep sex 

-2*fit_3$loglik[2]
# omition of race increases -2loglik so reduces model fit, therefore we want to keep sex 

-2*fit_4$loglik[2]
# omition of stage increases -2loglik so reduces model fit, therefore we want to keep sex 

# if something increases the -2 loglik of the model when omitted then we want to include this in our model 
# generally when we add variables to a model, if they decrease/reduce the -2loglik then they maximize the liklihood 

# test if there is a significant change in loglik between models 
anova(fit_2, fit_1) # P = 0.08
anova(fit_3, fit_1) # P = 0.02
anova(fit_4, fit_1) # P = 0.4
# exclusion of sex, race and stage all increase the -2loglik but only the increases in age and race are statistically 
# significant according to liklihood ratio test, there sex and race should definitely be included in the model 
# but stage doesn't necessairly need to be, but may be a good idea to keep it as we know it is predictive of survival
# in most data sets. can always compare models w/ + w/o stage in final analysis 

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

#------------------------------------------------------------------------------------------
# multivariate cox models 
#------------------------------------------------------------------------------------------

# fit model w/ mir_100 as categorial variable 
fit_cat <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                   covs_oral_2$mir_100
                 + covs_oral_2$age_at_initial_pathologic_diagnosis 
                 + covs_oral_2$gender 
                 + covs_oral_2$race_list
                 + covs_oral_2$tumor_stage)
# fit model w/ mir_100 as continuous variable 
fit_cont <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                    covs_oral_2$mir_100_continuous 
                  + covs_oral_2$age_at_initial_pathologic_diagnosis 
                  + covs_oral_2$gender 
                  + covs_oral_2$race_list
                  + covs_oral_2$tumor_stage)
# fit model w/ mir_100 as categorial variable w/o stage
fit_cat_sub <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                       covs_oral_2$mir_100
                     + covs_oral_2$age_at_initial_pathologic_diagnosis 
                     + covs_oral_2$gender 
                     + covs_oral_2$race_list)
#+ covs_oral_2$tumor_stage)
# fit model w/ mir_100 as continuous variable w/o stage
fit_cont_sub <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                        covs_oral_2$mir_100_continuous 
                      + covs_oral_2$age_at_initial_pathologic_diagnosis 
                      + covs_oral_2$gender 
                      + covs_oral_2$race_list)
#+ covs_oral_2$tumor_stage)

summary(fit_cat)
summary(fit_cont)
summary(fit_cat_sub)
summary(fit_cont_sub)
summary(fit_cat)$coef[1,5]
summary(fit_cat_sub)$coef[1,5]

summary(fit_cat)$coef[1,5]
exp(summary(fit_cat)$coef[1,1])
exp(confint.default(fit_cat))[1, ][1]
exp(confint.default(fit_cat))[1, ][2]

#save(expression, covs_oral_2, covs_3, file = "/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/5_TCGA_miRseq/Data_files/mirna_seq_surv_analysis.Rdata")
#save(covs_oral_2, file = "/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/5_TCGA_miRseq/Data_files/covs_oral_processed_surv_analysis.Rdata")

# build results table 
mir_100_results <- matrix(NA, nrow = 2, ncol = 8)
colnames(mir_100_results) <- c("Un_P", "Un_HR", "Un_95_CI_LB", "Un_95_CI_UB", 
                              "Ad_P", "Ad_HR", "Ad_95_CI_LB", "Ad_95_CI_UB")

# adjusted model 
mir_100_results[1, "Ad_P"] <- round(summary(fit_cat)$coef[1,5], digits = 8)
mir_100_results[1, "Ad_HR"] <- round(exp(summary(fit_cat)$coef[1,1]), digits = 8)
mir_100_results[1, "Ad_95_CI_LB"] <- round(exp(confint.default(fit_cat))[1, ][1], digits = 8)
mir_100_results[1, "Ad_95_CI_UB"] <- round(exp(confint.default(fit_cat))[1, ][2], digits = 8)

# fit unadjusted model w/ mir_100 as categorial variable 
fit_cat_unadj <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                         covs_oral_2$mir_100)

# unadjusted model
mir_100_results[1, "Un_P"] <- round(summary(fit_cat_unadj)$coef[1,5], digits = 8)
mir_100_results[1, "Un_HR"] <- round(exp(summary(fit_cat_unadj)$coef[1,1]), digits = 8)
mir_100_results[1, "Un_95_CI_LB"] <- round(exp(confint.default(fit_cat_unadj))[1, ][1], digits = 8)
mir_100_results[1, "Un_95_CI_UB"] <- round(exp(confint.default(fit_cat_unadj))[1, ][2], digits = 8)

# adjusted model 
mir_100_results[2, "Ad_P"] <- round(summary(fit_cont)$coef[1,5], digits = 8)
mir_100_results[2, "Ad_HR"] <- round(exp(summary(fit_cont)$coef[1,1]), digits = 8)
mir_100_results[2, "Ad_95_CI_LB"] <- round(exp(confint.default(fit_cont))[1, ][1], digits = 8)
mir_100_results[2, "Ad_95_CI_UB"] <- round(exp(confint.default(fit_cont))[1, ][2], digits = 8)

# fit unadjusted model w/ mir_100 as categorial variable 
fit_con_unadj <- coxph(Surv(covs_oral_2$time_to_event, covs_oral_2$vital_status_2) ~ 
                         covs_oral_2$mir_100_continuous)

# unadjusted model
mir_100_results[2, "Un_P"] <- round(summary(fit_con_unadj)$coef[1,5], digits = 8)
mir_100_results[2, "Un_HR"] <- round(exp(summary(fit_con_unadj)$coef[1,1]), digits = 8)
mir_100_results[2, "Un_95_CI_LB"] <- round(exp(confint.default(fit_con_unadj))[1, ][1], digits = 8)
mir_100_results[2, "Un_95_CI_UB"] <- round(exp(confint.default(fit_con_unadj))[1, ][2], digits = 8)

# CONCLUSIONS: 
# ~ mir_100 not associated as a continuous variable but is as a categorical variable 
# ~ mir_100 association as a catgorical variable not really affected by +/- stage in models 
setwd("/Users/Owen/GitHub/HNSCC_miRSNP/")
#write.table(mir_100_results, 'Validation/Analysis/7_TCGA_expression/Data_files/coxph_model_mir_100_exp.csv', 
            #sep=',' , col.names = T, row.names = F)
save(mir_100_results, file = "Validation/Analysis/7_TCGA_expression/Data_files/coxph_model_mir100_exp.Rdata")
# write CSV w/ average expression values accross all subjects for each miR of interest 
#write.table(avg_rpm_2, file = "Validation/Analysis/5_TCGA_miRseq/Data_files/miR_of_interest_TCGA_HNSCC_survival.csv", sep = ",", row.names = T, col.names = NA)


