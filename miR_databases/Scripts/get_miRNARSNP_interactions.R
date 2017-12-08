###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Curate predicted miRNA target sites predicted to be affected by HNSCC-associated SNPs
# Authors: Owen Wilkins 
# Link to miRNASNPv2 online database (http://bioinfo.life.hust.edu.cn/miRNASNP2/index.php)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")

#------------------------------------------------------------------------------------------
# read in data
#------------------------------------------------------------------------------------------

#utr_gain <- read.delim("HNSCC_miRSNP/MIRNASNP/miRNA_gain_by_SNPs_in_gene_3utr.txt", header = F, stringsAsFactors = FALSE)
#utr_loss <- read.delim("HNSCC_miRSNP/MIRNASNP/miRNA_loss_by_SNPs_in_gene_3utr.txt", header = F, stringsAsFactors = FALSE)
#save(utr_gain, file = "HNSCC_miRSNP/MIRNASNP/MIRNASNP_3utr_gain.txt")
#save(utr_loss, file = "HNSCC_miRSNP/MIRNASNP/MIRNASNP_3utr_loss.txt")

#------------------------------------------------------------------------------------------
# load data back in 
#------------------------------------------------------------------------------------------
load("miR_databases/Data/miRNASNP/MIRNASNP_3utr_gain.txt")
load("miR_databases/Data/miRNASNP/MIRNASNP_3utr_loss.txt")
head(utr_gain)
head(utr_loss)

#------------------------------------------------------------------------------------------
# extract entires for HNSCC-associated SNPs from target site gain and loss data sets 
#------------------------------------------------------------------------------------------

# write function to index dataframe based on SNP identifier 
find_interactions <- function(SNP_id, data_set){
  data_set_2 <- data_set[data_set$V4 == SNP_id,]
}

# apply function to identify miRNA:mRNA interactions 
snps <- c("rs41276823", "rs3831960", "rs16988668", 
          "rs75820821", 
          "rs56312243", "rs77506493", "rs45494801", 
          "rs56161233", "rs3777710", "rs12280753")
snp_gain <- lapply(snps, find_interactions, utr_gain)
snp_loss <- lapply(snps, find_interactions, utr_loss)

#------------------------------------------------------------------------------------------
# combine and save the results 
#------------------------------------------------------------------------------------------
comb_gain <- as.data.frame(rbind(snp_gain[[1]], snp_gain[[2]], snp_gain[[3]], 
                                 snp_gain[[4]], 
                                 snp_gain[[5]], snp_gain[[6]], snp_gain[[7]], 
                                 snp_gain[[8]], snp_gain[[9]], snp_gain[[10]])) 
comb_loss <- as.data.frame(rbind(snp_gain[[1]], snp_gain[[2]], snp_gain[[3]], 
                                 snp_gain[[4]], 
                                 snp_gain[[5]], snp_gain[[6]], snp_gain[[7]], 
                                 snp_gain[[8]], snp_gain[[9]], snp_gain[[10]])) 
# output results 
write.csv(comb_gain, file = "miR_databases/Survival/Output/survival_MIRNASNP_interactions_gain2.csv")
write.csv(comb_loss, file = "miR_databases/Survival/Output/survival_MIRNASNP_interactions_loss2.csv")

