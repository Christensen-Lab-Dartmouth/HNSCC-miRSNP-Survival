
###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TCHA HNSCC covariate data cleaning & model building for multivariable survival analyses including TCGA exp. data 
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen 1//Thesis/HNSCC_miRSNP/")

#------------------------------------------------------------------------------------------
# read data
#------------------------------------------------------------------------------------------
# covariate data 
covs <- readRDS("05_TCGA_expression_survival_analyses/Data_files/TCGA_HNSC_covariates_miRNA_expression_subset.rds")
covs$bcr_patient_uuid <- tolower(covs$bcr_patient_uuid)
# load subvject metadata 
metadata <- read.table("05_TCGA_expression_survival_analyses/Data_files/metadata.cart.2017-04-13T19-36-30.619378.txt", sep = "\t", header = T, stringsAsFactors = F)

#------------------------------------------------------------------------------------------
# clean covariates of interest & explore other variables 
#------------------------------------------------------------------------------------------

# combine death and follow uop variable to get time to event for alive and dead subjects 
covs$time_to_event <- NA
covs$time_to_event[is.na(covs$days_to_death)] <- covs$days_to_last_followup[is.na(covs$days_to_death)]
covs$time_to_event[!is.na(covs$days_to_death)] <- covs$days_to_death[!is.na(covs$days_to_death)]
covs$time_to_event

# consider age, sex, race and tumor stage variables 
covs$age_at_initial_pathologic_diagnosis
covs$gender
covs$gender <- factor(covs$gender, levels = c("MALE", "FEMALE"))
covs$race_list

# process race variable to 'white' 'non-white' 
covs$race_list[covs$race_list == "BLACK OR AFRICAN AMERICAN"] <- "NON-WHITE"
covs$race_list[covs$race_list == "ASIAN"] <- "NON-WHITE"
covs$race_list[covs$race_list == "BLACK OR AFRICAN AMERICAN"] <- "NON-WHITE"
covs$race_list[covs$race_list == ""] <- NA
covs$race_list

# remove those with NAs for race + re-level 
covs$race_list <- factor(covs$race_list, levels = c("WHITE", "NON-WHITE"))
covs_2 <- covs[!is.na(covs$race_list),]

# split complicated stage variable up and extract stage number data 
x <- strsplit(covs_2$stage_event, " ")
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
covs_2$tumor_stage <- tumor_stage
covs_2$tumor_stage <- factor(covs_2$tumor_stage, levels = c("I+II", "III+IV"))

saveRDS(covs_2, file = "05_TCGA_expression_survival_analyses/Data_files/TCGA_HNSC_covariates_processed.rds")

