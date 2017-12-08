
###########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Process TCGA-HNSCC subject metadata file from JSON files to txt file format 
# Authors: Owen Wilkins 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################
rm(list = ls())
library(rjson)
library(plyr)
setwd("/Users/Owen/Thesis/HNSCC_miRSNP/")

# set file paths for in and out files 
myinf1 = "TCGA_RNAseq_analysis/Data_files/TCGA_HNSC_metadata_2017-04-07T22-01-38.673662.json"
myoutf1 = "TCGA_RNAseq_analysis/Data_files/TCGA_HNSC_metadata_2017-04-07T22-01-38.673662.txt"

# convert JSON object into R object 
mydata = fromJSON(file=myinf1)

# There are some list values that are NULL; replace with NA
myfun = function(arg) {
  arg[is.null(arg)] = NA
  return(arg)
}
mydata = rapply(mydata, f=myfun, how="list")

# Extract values from list; output will have names that will indicate the list level they were originally part of
res = unlist(mydata[[1]], recursive=T)
mynam = names(res)

# One patient may have multiple sets of records. Therefore, you must remove duplicates.
res = res[!duplicated(mynam)] 
res = t(res)
res = as.data.frame(res)

for (i in 2:length(mydata)) {
  cat("\r", i)
  tmp = unlist(mydata[[i]], recursive=T)
  mynam = names(tmp)
  tmp = tmp[!duplicated(mynam)]
  tmp = t(tmp)
  tmp = as.data.frame(tmp)
  res = rbind.fill(res, tmp)
}

# write out processed metadata file 
write.table(res, myoutf1, sep="\t", row.names=F, quote=F)



