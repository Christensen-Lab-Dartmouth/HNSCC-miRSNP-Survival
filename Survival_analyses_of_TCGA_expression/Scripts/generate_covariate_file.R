###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate a single covariate file for all the TCGA-HNSCC clinical data that was downloaded from GDC 
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
library(XML)
setwd("TCGA_RNAseq_analysis/Clinical_data/")

# move files into one folder
clinical_folders = list.files()
#files = clinical_folders[grep("gdc", clinical_folders)]
#clinical_folders = clinical_folders[!(clinical_folders %in% files)]
if(!dir.exists("TCGA-HNSC_READ_clinical"))dir.create("TCGA-HNSC_READ_clinical")

# Process XML files 
setwd('TCGA-HNSC_READ_clinical')
clin_files <- list.files()
df <- xmlToDataFrame(clin_files[1])
tag <- is.na(df[1,])
df[1,tag] <- df[2,tag]
df <- df[1,]
row.names(df) <- df$bcr_patient_barcode
gdc_clin <- df

for (i in 2:length(clin_files)){
  df <- xmlToDataFrame(clin_files[i])
  tag <- is.na(df[1,])
  df[1,tag] <- df[2,tag]
  df <- df[1,]
  row.names(df) <- df$bcr_patient_barcode
  gdc_clin <- rbind(gdc_clin,df)
}

# write out csv files 
write.csv(gdc_clin, file = 'TCGA_RNAseq_analysis/TCGA_HNSC_covariates.csv')

