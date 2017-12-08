###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get miRNAseq isoform specific expression data for miRNAs of interest from TCGA adjacent normal tissue samples from HNSCC subjects 
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())

#------------------------------------------------------------------------------------------------------------------------------
# Load + prepare data 
#------------------------------------------------------------------------------------------------------------------------------

setwd("TCGA_miRseq_analysis/Data_files/TCGA_HNSC_miRNA-seq_isoforms_normal_tissue/")
dir()

# get all file name in directory 
filelist = list.files(path = dir(),pattern = "*.*.txt")

# remove subjects with neoadjuvant tx or were identified as recurrent tumors
# these sujects are indcated by annotation files 
indicies_to_drop <- grep("annotations.txt", filelist)
filelist_2 <- filelist[-indicies_to_drop]

# run a loop over elements of file list and make a new vector to add file path to them 
paths <- c()
for(i in 1:length(filelist_2)){
  paths[i] <- paste(dir()[i], "/", filelist_2[i], sep = "")
  paths
}

# read in all files 
datalist = lapply(paths, function(x)read.table(x, header=T, stringsAsFactors = F)) 

# bind all rows of these data frames together 
mirseq = do.call("rbind", datalist)
summary(mirseq$reads_per_million_miRNA_mapped)
##### note: that there are no 0.0000 values for any entries in this isofrm specific data 
##### (which is not the case with the averaged data)

#------------------------------------------------------------------------------------------------------------------------------
# extract and save average reads per million (rpm) for each miR of interest 
#------------------------------------------------------------------------------------------------------------------------------
length(unique(mirseq$miRNA_region)) # check how many unique miRNA entries exist using miRBase identifier

# calculate mean expression values for all miRs 
calc_mean_exp <- function(data_set){
  mir_mean_exp <- data.frame() # estalish data frame to fill 
  unique_miRs <- unique(data_set$miRNA_region) # get vector of unique miRNA miRbase identifiers
  for(i in 1:length(unique_miRs)){
    # loop over each unique miRBase identifier, create data frame of expression values for these entries
    # average expression and add additional info on each entry 
    exp_data <- data_set[which(data_set$miRNA_region == unique_miRs[i]),]
    mir_mean_exp[i, 1] <- exp_data$miRNA_ID[1]
    mir_mean_exp[i, 2] <- mean(exp_data$read_count)
    mir_mean_exp[i, 3] <- log2(mean(exp_data$read_count))
    mir_mean_exp[i, 4] <- mean(exp_data$reads_per_million_miRNA_mapped)
    mir_mean_exp[i, 5] <- log2(mean(exp_data$reads_per_million_miRNA_mapped))
    mir_mean_exp[i, 6] <- names(table(exp_data$miRNA_region))[1]
  }
  # remove any entries which are stem loop, precursor or unannotated 
  mir_mean_exp <- mir_mean_exp[-which(mir_mean_exp[,6] == "precursor" | mir_mean_exp[,6] == "stemloop" | 
                                        mir_mean_exp[,6] == "unannotated"),]
  colnames(mir_mean_exp) <- c("miR", "mean_read_ct", "log2_mean_read_ct", "mean_rpm", "log2_mean_rpm", "MIMAT")
  mir_mean_exp
}

# apply function to calculate 
avg_rpm <- calc_mean_exp(mirseq)

# check distributions and summaries 
hist(avg_rpm$mean_rpm)
hist(log2(avg_rpm$mean_rpm))
summary(avg_rpm$mean_rpm)
summary(log2(avg_rpm$mean_rpm))

# look at quantile distribution 
quantile(avg_rpm$mean_rpm, probs = c(0, 0.25, 0.5, 0.75, 1)) 
quantile(log2(avg_rpm$mean_rpm), probs = c(0, 0.25, 0.5, 0.75, 1)) 

# calculate percentile rank using 
perc.rank <- ecdf(log2(avg_rpm$mean_rpm))
avg_rpm$percentile <- sapply(log2(avg_rpm$mean_rpm), perc.rank)

#------------------------------------------------------------------------------------------------------------------------------
# subset expression data for miRNAs of interest 
# (note: indexing performed using miRBase identifiers, so that correct miRNA strand is identified)
#------------------------------------------------------------------------------------------------------------------------------

# read in list of miRs of interest with mirbase identifiers 
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")
#mimats_of_int <- read.csv("/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/5_TCGA_miRseq/mirna_mirbase_accession_links_survival_miRNA_name_only.csv", 
#header = TRUE, stringsAsFactors = FALSE)
#mimats_of_int <- mimats_of_int[-which(duplicated(mimats_of_int)),]
#write.csv(mimats_of_int, "/Users/Owen/GitHub/HNSCC_miRSNP/Validation/Analysis/5_TCGA_miRseq/mirna_mirbase_accession_links_survival_unique_miRs.csv", row.names = F)
mimats_of_int <- read.csv("TCGA_miRseq_analysis/Survival/Data_files/mirna_mirbase_accession_links_survival_unique_miRs_w_miRBase_IDs.csv", 
                          header = TRUE, stringsAsFactors = FALSE) 

# removes repeated miRbase identifers (MIMATS)
mimats_of_int_2 <- mimats_of_int[!duplicated(mimats_of_int$MIMAT),]
table(duplicated(mimats_of_int_2$mirNA_name))

# remove "mature," from the character string for mirbase identifiers in mirseq data 
for(i in 1:length(avg_rpm$MIMAT)){
  avg_rpm$MIMAT_2[i] <- strsplit(avg_rpm$MIMAT, ",")[[i]][2]
}

# how many MIMATs of interest are found in the averaged isform data 
table(avg_rpm$MIMAT_2 %in% unique(mimats_of_int_2$MIMAT))

# index data for only these miRs
avg_rpm_2 <- avg_rpm[avg_rpm$MIMAT_2 %in% unique(mimats_of_int_2$MIMAT),] 

# Note that some MIMATs were not found, so look at them here to check them out 
mimats_not_found <- mimats_of_int_2[!unique(mimats_of_int_2$MIMAT) %in% avg_rpm$MIMAT_2,]
length(unique(mimats_not_found$MIMAT))

# write out averaged exprsssion data  
write.csv(avg_rpm_2, file = "TCGA_miRseq_analysis/Survival/Data_files/TCGA_HNSCC_miRNA_exp_normal_tissue_isoforms_survival.csv")


