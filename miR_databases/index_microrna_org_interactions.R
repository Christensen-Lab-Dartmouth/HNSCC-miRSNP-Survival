###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify predicted microRNA target sites from microrna.org that overlap with overall survival associated SNPs in HNSCC
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")

#------------------------------------------------------------------------------------------
# read in data sets 
#------------------------------------------------------------------------------------------
# load downloads for conserved sites with good mirSVR scores and the same for non-conserved sites 
#mir_0_0_scores <- read.delim("miR_databases/Data/human_predictions_0_0_aug2010.txt", header = T, stringsAsFactors = FALSE)
#mir_0_C_scores <- read.delim("miR_databases/Data/human_predictions_0_C_aug2010.txt", header = T, stringsAsFactors = FALSE)
#mir_S_0_scores <- read.delim("miR_databases/Data/human_predictions_S_0_aug2010(1).txt", header = T, stringsAsFactors = FALSE)
#mir_S_C_scores <- read.delim("miR_databases/Data/human_predictions_S_C_aug2010(1).txt", header = T, stringsAsFactors = FALSE)

# join into list 
#all_scores <- rbind(mir_0_0_scores, mir_0_C_scores, mir_S_0_scores, mir_S_C_scores)

# save output 
#save(all_scores, file = "miR_databases/Data/all_scores_list.Rdata")

#--------------------------------------------------------------------------------------------------------------------------
# load saved data
#--------------------------------------------------------------------------------------------------------------------------
load("miR_databases/Data/microrna.org/all_scores_list.Rdata")

# save a histogram of mirsvr score distribution 
dpi = 300
png("Validation/Analysis/4_miR_databases/Data/mirsvr_distribution.png", width = dpi*7, height = dpi*8, res = dpi)
hist((all_scores$mirsvr_score*-1), xlab = "Reciprocal mirSVR scores", ylab = "No. of predicted target sites")
abline(v =0.1, lty = 2, lwd = 2, col = "red")
dev.off()

# calculate percentile ranks for mirsvr scores and add to all_scores
perc.rank <- ecdf(all_scores$mirsvr_score)
all_scores$percentile <- sapply(all_scores$mirsvr_score, perc.rank)
head(all_scores)

#--------------------------------------------------------------------------------------------------------------------------
# extract genomic coordinates 
#--------------------------------------------------------------------------------------------------------------------------

# check gene identifiers for all genes containing HNSCC associated SNPs are in data set 
genes <- c("RGS3", "ZSWIM5", "ZSCAN22", "CCDC97", "RPL28", "ZNF766", "SH3BP4", "FAM8A1", "BUD13")
genes_indicies <- sapply(genes, function(x) which(all_scores$gene_symbol==x))
str(genes_indicies)

# extract all predictions for each of these genes 
all_scores_genes <- lapply(genes_indicies, function(x) all_scores[x,])
str(all_scores_genes)

# write a function to extract genomic coordinates
extract_genome_coords <- function(data_set){
  # make vectors of NAs to fill for CHR and miR binding sites 
  # note some transcripts have more than one predicted binding site for a miR, therefore need 2 vectors for these 
  data_set$CHR <- NA
  data_set$site_start_1 <- NA
  data_set$site_end_1 <- NA
  data_set$site_start_2 <- NA
  data_set$site_end_2 <- NA
  # loop over rows in data set and extract genome coordinates 
  for(i in 1:nrow(data_set)){
    # extract chromosome and fill CHR 
    data_set$CHR[i] <- as.numeric(strsplit(data_set$genome_coordinates[i], ":")[[1]][2])
    # split and extract position from genome coordinates 
    coord_split <- strsplit(data_set$genome_coordinates[i], ":")[[1]][3]
    coord_split_2 <- strsplit(coord_split, ",")[[1]]
    # if there is more than 1 predicted site, fill both vectors for miR binding sites
    if(length(coord_split_2) > 1){
      coord_split_3 <- strsplit(coord_split_2, "-")
      data_set$site_start_1[i] <- as.numeric(coord_split_3[[1]][1])
      data_set$site_end_1[i] <- as.numeric(coord_split_3[[1]][2])
      data_set$site_start_2[i] <- as.numeric(coord_split_3[[2]][1])
      data_set$site_end_2[i] <- as.numeric(coord_split_3[[2]][2])
      # otherwise, just fill the 1st set of vectors 
    } else {
      data_set$site_start_1[i] <- as.numeric(strsplit(coord_split_2, "-")[[1]][1])
      data_set$site_end_1[i] <- as.numeric(strsplit(coord_split_2, "-")[[1]][2])
    }
  }
  data_set
}

# apply function to extract all genomic coordinates  
all_scores_genes_2 <- lapply(all_scores_genes, extract_genome_coords)
str(all_scores_genes_2)

#--------------------------------------------------------------------------------------------------------------------------
# index miR data for only interactions (rows) intersecting HNSCC-associated SNPs
#--------------------------------------------------------------------------------------------------------------------------

# write function to identify miRNA-mRNA predicted interactions overlapping SNP loci 
find_interactions <- function(data_set, chromosome, SNP_loci){
  # create a new column that can be used as an indicator variable for if an interactions meets the requirements
  data_set$rows_test_1 <- rep(0, dim(data_set)[1])
  data_set$rows_test_2 <- rep(0, dim(data_set)[1])
  # loop over all rows in data_set and check if the SNP_loci falls into the genomic range of the miR 
  for(i in 1:dim(data_set)[1]){
    # generate a vector of the genomic range spanned by the miR 
    if(is.na(data_set$site_start_2[i])){
      ranges <- seq(data_set$site_start_1[i], data_set$site_end_1[i])
      # if SNP loci is encompassed by miR loci, put a 1 in rows_test 
      if(SNP_loci %in% ranges){
        data_set$rows_test_1[i] <- 1
      }
    } else{
      ranges_1 <- seq(data_set$site_start_1[i], data_set$site_end_1[i])
      ranges_2 <- seq(data_set$site_start_2[i], data_set$site_end_2[i])
      if(SNP_loci %in% ranges_1){
        data_set$rows_test_1[i] <- 1
      }
      if(SNP_loci %in% ranges_2){
        data_set$rows_test_2[i] <- 1
      }
    }
    # index data frame for only columns with a 1
    data_set_2 <- data_set[data_set$rows_test_1 == 1 | data_set$rows_test_2 == 1,]
  }
  # order data frame by mirSVR score 
  data_set_3 <- data_set_2[order(data_set_2$mirsvr_score),]
  # return our new dataset 
  data_set_3
}

# apply function to locate interactaions encompassing coordinates of SNPs of interest 
# overall HNSCC
rs41276823_ints <- find_interactions(all_scores_genes_2$RGS3, 9, 116359270)
rs3831960_ints <- find_interactions(all_scores_genes_2$ZSWIM5, 1, 45483218)
rs16988668_ints <- find_interactions(all_scores_genes_2$ZSCAN22, 19, 58850756)
# oral cavity cancer
rs75820821_ints <- find_interactions(all_scores_genes_2$CCDC97, 19, 41829457)
# pharyngeal cancer
rs56312243_ints <- find_interactions(all_scores_genes_2$RPL28, 19, 55899602)
rs77506493_ints <- find_interactions(all_scores_genes_2$ZNF766, 19, 52795158)
rs45494801_ints <- find_interactions(all_scores_genes_2$RPL28, 19, 55902084)
# laryngeal cancer 
rs56161233_ints <- find_interactions(all_scores_genes_2$SH3BP4, 2, 235964133)
rs3777710_ints <- find_interactions(all_scores_genes_2$FAM8A1, 6, 17611089)
rs12280753_ints <- find_interactions(all_scores_genes_2$BUD13, 11, 116613660)

# index interaction data for each SNP for only mirsvr score >-0.1
index_mirsvr <- function(data_set){
  data_set_2 <- data_set[as.numeric(data_set$mirsvr_score) < -0.1,]
  data_set_2
}

# collate into list 
interactions <- list(rs41276823_ints, rs3831960_ints, rs16988668_ints, 
                     rs75820821_ints, 
                     rs56312243_ints, rs77506493_ints, rs45494801_ints, 
                     rs56161233_ints, rs3777710_ints, rs12280753_ints)

# apply indexing for mirsvr scores over list 
interactions_sub <- lapply(interactions, index_mirsvr)
str(interactions_sub)

# bind into 1 data frame 

interactions_sub_2 <- as.data.frame(rbind(interactions_sub[[1]], interactions_sub[[2]], interactions_sub[[3]], 
                            interactions_sub[[4]], 
                            interactions_sub[[5]], interactions_sub[[6]], interactions_sub[[7]],
                            interactions_sub[[8]], interactions_sub[[9]], interactions_sub[[10]]))

# save output tables 
write.table(interactions_sub_2, file = "miR_databases/Survival/Output/survival_microrna_org_high_confidence_interations.csv", sep = ",", row.names = T, col.names = NA)
