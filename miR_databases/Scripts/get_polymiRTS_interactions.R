###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract data for predicted miRNA target sites at HNSCC-associated SNPs 
# Authors: Owen Wilkins 
# Link: (http://compbio.uthsc.edu/miRSNP/)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/GitHub/HNSCC_miRSNP/")

#------------------------------------------------------------------------------------------
# read in data
#------------------------------------------------------------------------------------------
targets <- read.delim("miR_databases/Data/PolymiRTS/target_miRSNP_human.txt", header = T, stringsAsFactors = FALSE)

#------------------------------------------------------------------------------------------
# extract data for HNSCC-associated SNPs 
#------------------------------------------------------------------------------------------

# write function to index targets dataframe based on SNP ID 
find_interactions <- function(SNP){
  interactions <- list()
  # loop over SNP and CHR/POS to curate info on SNPs, and fill interactions list 
  for(i in 1:length(SNP)){
      interactions_temp <- as.data.frame(targets[targets$SNPID == SNP[i],])
      # add conditional to fill interaction list w/ NULL when no SNPs are found 
      if(dim(interactions_temp)[1] != 0){
        interactions[[i]] <- interactions_temp
        names(interactions)[[i]] <- SNP[i]
      } else {interactions[[i]] <- NULL}
  }
  interactions
}

# apply function to identify miR:mRNA interactions 
targets_ints <- find_interactions(c("rs41276823", "rs3831960", "rs16988668", 
                                    "rs75820821",
                                    "rs56312243", "rs77506493", "rs45494801",
                                    "rs56161233", "rs3777710", "rs12280753"))
# look at results 
str(targets_ints)
targets_ints

#------------------------------------------------------------------------------------------
# clean data in usable format 
#------------------------------------------------------------------------------------------

# write function to process data from each SNP of interest into a usable format
process <- function(data_set) {
  # make results tables to store data in 
  results_A1 <- data.frame()
  results_A2 <- data.frame()
  # split up vectors of interest 
  split_vectors <- function(vector_to_split) strsplit(vector_to_split, "|", fixed = TRUE)
  A1_miRs <- sapply(data_set[, "Allele1miR"], split_vectors)
  A1_fun <- sapply(data_set[, "Allele1Func"], split_vectors)
  A1_cs_diff <- sapply(data_set[, "cs_diff_Allele1Site"], split_vectors)
  # for the other allele too 
  A2_miRs <- sapply(data_set[, "Allele2miR"], split_vectors)
  A2_fun <- sapply(data_set[, "Allele2Func"], split_vectors)
  A2_cs_diff <- sapply(data_set[, "cs_diff_Allele2Site"], split_vectors)
  # iterate over the new items and extract related data 
  for(i in 1:length(A1_miRs[[1]])){
    results_A1[i, 1] <- data_set$Allele1
    results_A1[i, 2] <- A1_miRs[[1]][i]
    results_A1[i, 3] <- A1_fun[[1]][i]
    results_A1[i, 4] <- A1_cs_diff[[1]][i]
  }
  for(i in 1:length(A2_miRs[[1]])){
    results_A2[i, 1] <- data_set$Allele2
    results_A2[i, 2] <- A2_miRs[[1]][i]
    results_A2[i, 3] <- A2_fun[[1]][i]
    results_A2[i, 4] <- A2_cs_diff[[1]][i]
  }
  results <- rbind(results_A1, results_A2)
  colnames(results) <- c("Allele", "miR", "function_class", "cs_score_diff")
  results
}
# apply this for each SNP 
targets_ints_2 <- lapply(targets_ints, process)

# combine them together into one data frame 
comb <- as.data.frame(rbind(targets_ints_2[[1]], targets_ints_2[[2]], targets_ints_2[[3]], 
                            targets_ints_2[[4]], 
                            targets_ints_2[[5]], targets_ints_2[[6]], targets_ints_2[[7]], 
                            targets_ints_2[[8]], targets_ints_2[[9]]))
rm(targets_ints, targets_ints_2)

#------------------------------------------------------------------------------------------
# calculate and add in context+ score percentiles 
#------------------------------------------------------------------------------------------

# get all context+ scores
split_vectors <- function(vector_to_split) strsplit(vector_to_split, "|", fixed = TRUE)
A1_cs_diffs <- na.omit(as.numeric(unlist(sapply(targets[, "cs_diff_Allele1Site"], split_vectors))))
A2_cs_diffs <- na.omit(as.numeric(unlist(sapply(targets[, "cs_diff_Allele2Site"], split_vectors))))
cs_diffs <- c(A1_cs_diffs, A2_cs_diffs)

# look at distribution 
hist(cs_diffs)
quantile(cs_diffs, probs = c(0, 0.25, 0.5, 0.75, 1)) 

# calculate percentiles for cs_diff scores 
perc.rank <- ecdf(cs_diffs)
percentile <- sapply(cs_diffs, perc.rank)

# calculate and add context+ score percentiles to interaction summary table 
str(comb)
comb$cs_score_diff <- as.numeric(comb$cs_score_diff)
comb$cs_score_diff_perc <- sapply(comb$cs_score_diff, perc.rank)

write.csv(comb, file = "miR_databases/Survival/Output/survival_polymirts_interactions.csv")
