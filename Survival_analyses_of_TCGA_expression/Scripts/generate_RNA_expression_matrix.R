###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in miRNAseq data from TCGA HNSCC cases and generate expression matrix
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")

#------------------------------------------------------------------------------------------
# read in data
#------------------------------------------------------------------------------------------

# read in covaraite file 
covs <- readRDS("TCGA_expression_survival_analyses/Data_files/TCGA_HNSC_covariates_miRNA_expression_subset.rds")
covs$bcr_patient_uuid <- tolower(covs$bcr_patient_uuid)

# load in metadata 
metadata <- read.table("TCGA_expression_survival_analyses/Data_files/TCGA_HNSC_metadata_2017-04-07T22-01-38.673662.txt", sep = "\t", header = T, stringsAsFactors = F)

# RNA_seq data
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/TCGA_expression_survival_analyses/RNA-seq/")
dir("/Users/Owen/Thesis/HNSCC_miRSNP/TCGA_expression_survival_analyses/RNA-seq/")

#------------------------------------------------------------------------------------------
# read in and process RNA-seq data for analysis 
#------------------------------------------------------------------------------------------

# get all file name in directory 
filelist = list.files(path = dir(),pattern = ".txt", full.names = TRUE)
filelist

# note: there are "annotation.txt" files associated with 36 subjects that indicate if case was found to be a recurrence
# after submission, if the patient had a simulatneous cancer of another kind etc. 
# these cases should be removed 
indicies_to_drop <- grep("annotations.txt", filelist)
filelist_2 <- filelist[-indicies_to_drop]

# read in all files 
datalist = lapply(filelist_2, function(x) read.table(x, header=F, stringsAsFactors = F)) 

# split the IDs from filelist_2 and match them to the metadata 
get_ids <- function(vector_to_split) strsplit(vector_to_split, "/")[[1]][2]
filelist_3 <- sapply(filelist_2, get_ids)
names(datalist) <- filelist_3
names(datalist)[1:5]
names(filelist_3) <- NULL

# check that all filelist items are in metadata 
table(filelist_3 %in% metadata$file_name)
# 2 is missing as filelist has 502 eklements and metadata data has 500

filelist_4 <- filelist_3[match(metadata$file_name, filelist_3)] # index filelist for subjects w/ metadata 
filelist_5 <- filelist_4[order(filelist_4)] # order metadata and filelist the same way  
metadata_2 <- metadata[order(metadata$file_name),]
all(metadata_2$file_name == filelist_5) # check they are identical 

datalist_2 <- datalist[match(filelist_5, names(datalist))] # index datalist for elements w/ matching IDs in filelist_5
datalist_3 <- datalist_2[order(names(datalist_2))] # order datalist_2 same way as filelist_5
all(filelist_5 == names(datalist_3)) # check they are identical 

names(datalist_3)[1:5]
filelist_5[1:5]
metadata_2$file_name[1:5]

names(datalist_3) <- metadata_2$cases.case_id # give datalist patient_uuid as names (instead of file uuid)
datalist_4 <- datalist_3[order(names(datalist_3))] # order based on this new variable 

#------------------------------------------------------------------------------------------
# add patient UUID to RNA_seq data (to link expression and clinical data)
#------------------------------------------------------------------------------------------

# check to see if all subjects in covs are in the datalist 
table(tolower(covs$bcr_patient_uuid) %in% names(datalist_4))
# there are 5 missing, so remove those 5 when you index covariates for subjects with metadata 

# index covs for those in metadata 
covs_2 <- covs[na.omit(match(names(datalist_4), tolower(covs$bcr_patient_uuid))),] 
# index datalist for subjects in covariate data so they match 
datalist_5 <- datalist_4[na.omit(match(tolower(covs_2$bcr_patient_uuid), (names(datalist_4))))]
# order covariate data and datalist in same way 
covs_3 <- covs_2[order(tolower(covs_2$bcr_patient_uuid)),] 
datalist_6 <- datalist_5[order(names(datalist_5))]
# check datalist and covs_2 are in same order with same identifiers 
all(tolower(covs_3$bcr_patient_uuid) == names(datalist_6)) 

# make matrix to hold expression data 
expression <- matrix(NA, nrow = length(covs_3$bcr_patient_uuid), ncol = length(datalist_4[[1]]$V1)) 
colnames(expression) <- datalist_4[[1]]$V1
rownames(expression) <- covs_2$bcr_patient_uuid
for(i in 1:length(datalist_6)){
  expression[i,] <- datalist_6[[i]]$V2
} # fill matrix w/ expression data 
expression[1:10, 1:5]

# clean up workspace
rm(filelist_5, filelist_4, filelist_3, filelist_2, filelist, covs_2, covs, metadata, metadata_2, 
   datalist, datalist_2, datalist_3, datalist_4, indicies_to_drop)

#------------------------------------------------------------------------------------------
# subset subjects by expression of genes of interest + run Kaplan Meier analysis 
#------------------------------------------------------------------------------------------

# check indicies in expression matrix for genes of interest to see if present
which(colnames(expression) == "ENSG00000130147.15") # SH3BP4 # BOG25 # TPP # EHB10
which(colnames(expression) == "ENSG00000259571.1") # BLID

# check to see if all gene identifiers are ensembl gene IDs (so you know that is the only type of 
# identifier SH3BP4 might be found under in this data set)
x <- strsplit(colnames(expression), "")
y <- c()
for(i in 1:length(x)){
  y[i] <- c(paste(x[[i]][1], x[[i]][2], x[[i]][[3]], x[[i]][[4]]))
}
y
table(y == c("E N S G"))
dim(expression)
# complete match, so they do all seem to be ensembl gene IDs 

# save expression matrix for RNA-seq data 
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")
saveRDS(expression, file = "TCGA_expression_survival_analyses/Data_files/TCGA_HNSCC_tumor_RNAseq_expression_matrix.rds")
saveRDS(covs_3, file = "TCGA_expression_survival_analyses/Data_files/TCGA_oral_cancer_covs_cleaned_w_BLID_expression.rds")

rm(datalist_5, datalist_6, i, x, y, get_ids)

#------------------------------------------------------------------------------------------
# clean and pre-process covariate data from subset of subjects w/ RNA-seq data in the same way 
# as done for subjects with miRNAseq data in 'generate_miRNA_expression_matrix.R'
#------------------------------------------------------------------------------------------
covs <- covs_3
rm(covs_3)

# confirm subject order still matches 
all(covs$bcr_patient_uuid == rownames(expression))

# add expression data to covs 
covs$blid_continuous <- log2(expression[,which(colnames(expression) == "ENSG00000259571.1")])

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

covs_oral_3 <- covs_oral_2[-which(covs_oral_2$blid_continuous==-Inf),]

saveRDS(covs_oral_3, file = "TCGA_expression_survival_analyses/Data_files/TCGA_oral_cancer_covs_cleaned_w_BLID_expression.rds")
