###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify predicted microRNA target sites from MirSNP database that overlap with overall survival associated SNPs in HNSCC
# Authors: Owen Wilkins 
# Link: # (http://bioinfo.bjmu.edu.cn/mirsnp/search/)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")

#-----------------------------------------------------------------------------------------------------------------
# read in data 
#-----------------------------------------------------------------------------------------------------------------

#MirSNP <- read.delim("miR_databases/Data/MirSNPInTarget.txt", header = F, stringsAsFactors = FALSE)
#save(MirSNP, file = "miR_databases/Data/MirSNPInTarget.Rdata")
load("miR_databases/Data/mirSNP/MirSNPInTarget.Rdata")

#-----------------------------------------------------------------------------------------------------------------
# extract interactions for HNSCC-associated SNPs  
#-----------------------------------------------------------------------------------------------------------------

# provide column names using info from website 
colnames(MirSNP) <- c("rsid", "RefGene", "miRNA", "mRNA", "alleles", "mirSVR", "allele_changes", "impact_A1_A2", "A1", "A1_binding_score",
                      "FE_A1", "start_pos_A1", "end_pos_A1", "conservation_A1", "A2", "A2_binding_score",
                      "FE_A2", "start_pos_A2", "end_pos_A2", "conservation_A2")
head(MirSNP)

# specify SNPs of interest 
snps <- c("41276823", "3831960", "16988668", 
          "75820821",
          "56312243", "77506493", "45494801",
          "56161233", "3777710", "12280753")

# index to pull out SNPs of interest 
interactions <- lapply(snps, function(x) MirSNP[MirSNP$rsid == x,])

# combine together into 1 data frame 
comb <- as.data.frame(rbind(interactions[[1]], interactions[[2]], interactions[[3]], 
                            interactions[[4]], 
                            interactions[[5]], interactions[[6]], interactions[[7]], 
                            interactions[[8]], interactions[[9]], interactions[[10]]))

# write results to file 
write.table(comb, file = "miR_databases/Survival/Output/MirSNP_interactions_survival.csv", sep = ",", row.names = T, col.names = NA)
