###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in miRNAseq data from TCGA HNSCC cases and generate expression matrix
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")
library(survival)
library(survminer)

#------------------------------------------------------------------------------------------
# read in and clean up data
#------------------------------------------------------------------------------------------
# covariate data 
covs <- read.csv('TCGA_expression_survival_analyses/Data_files/TCGA_HNSC_covariates.csv', header = T, stringsAsFactors = F)
covs$bcr_patient_uuid <- tolower(covs$bcr_patient_uuid)

# load subvject metadata 
metadata <- read.table("TCGA_expression_survival_analyses/Data_files/metadata.cart.2017-04-13T19-36-30.619378.txt", sep = "\t", header = T, stringsAsFactors = F)

# set directory for miRseq data (not isoform specific)
setwd("TCGA_miRseq_analysis/Survival/Data_files/TCGA_HNSC_miRNA-seq/")
dir("TCGA_miRseq_analysis/Survival/Data_files/TCGA_HNSC_miRNA-seq/")

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
table(metadata$file_id %in% filelist_3)
# all TRUE 

# check if the cases.case_id variable in metadata is the patient uuid 
table(metadata$cases.case_id %in% covs$bcr_patient_uuid)
# all TRUE 

#covs$bcr_patient_uuid[covs$bcr_patient_barcode=="TCGA-4P-AA8J"]
#metadata$cases.case_id[metadata$cases.submitter_id=="TCGA-4P-AA8J"]
#covs$bcr_patient_uuid[covs$bcr_patient_barcode=="TCGA-F7-A50G"]
#metadata$cases.case_id[metadata$cases.submitter_id=="TCGA-F7-A50G"]
#covs$bcr_patient_uuid[covs$bcr_patient_barcode=="TCGA-4P-AA8J"]
#metadata$cases.case_id[metadata$cases.submitter_id=="TCGA-4P-AA8J"]
#covs$bcr_patient_uuid[covs$bcr_patient_barcode=="TCGA-UF-A7JF"]
#metadata$cases.case_id[metadata$cases.submitter_id=="TCGA-UF-A7JF"]

#filelist_4 <- filelist_3[match(metadata_2$file_name, filelist_3)] # index filelist for subjects w/ metadata 
filelist_4 <- filelist_3[order(filelist_3)] # order metadata and filelist the same way  
metadata_2 <- metadata[order(metadata$file_id),]
head(filelist_4) ; head(metadata_2$file_id)
# check they are identical
all(metadata_2$file_id == filelist_4) 

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

# save 
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")
saveRDS(covs_3, file = "TCGA_expression_survival_analyses/Data_files/TCGA_HNSC_covariates_miRNA_expression_subset.rds")
saveRDS(expression, file = "TCGA_expression_survival_analyses/Data_files/TCGA_HNSCC_tumor_miRNA_expression_matrix.rds")

