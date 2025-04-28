require(data.table)
library(tidyr)
require(tidyverse)

#load file
BrS <- fread("../data/Barc_2022_BrS_unclean.txt")

#clean file:

#delete columns that are not necessary
print(paste("BrS table delete columns:", nrow(BrS)))
BrS <- BrS %>% select(-chr, -pos, -n_cohorts, -n_cases, -n_controls, -direction, -hetISq, -hetChiSq, -hetDf, -hetPVal, -eaf_cases, -eaf_controls)
print(paste("BrS table delete columns(no difference!):", nrow(BrS)))

#rename columns
print(paste("BrS table rename:", nrow(BrS)))
colnames(BrS) <- c("SNP", "A1", "A2", "beta", "SE", "p_value", "N", "EAF") 
print(paste("BrS table rename(no difference!):", nrow(BrS)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(BrS$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  BrS$A1_multiple <- nchar(BrS$A1) > 1
  # delete those SNPs
  BrS <- BrS[BrS$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("BrS table delete A1_multiple:", nrow(BrS)))
  BrS <- BrS %>% select(-A1_multiple)
  print(paste("BrS table delete A1_multiple(no difference!):", nrow(BrS)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(BrS$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  BrS$A2_multiple <- nchar(BrS$A2) > 1
  # delete those SNPs
  BrS <- BrS[BrS$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("BrS table delete A2_multiple:", nrow(BrS)))
  BrS <- BrS %>% select(-A2_multiple)
  print(paste("BrS table delete A2_multiple(no difference!):", nrow(BrS)))
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- BrS$SNP == "" | !grepl("^r", BrS$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the once that are empty
  BrS <- BrS[!not_r_or_empty, ]
}

# Delete SNPs that are duplicates
print(paste("bevore delete duplicates:", nrow(BrS)))
BrS <- BrS[!duplicated(BrS$SNP), ]
print(paste("after delete duplicates:", nrow(BrS)))

#make MAF
BrS$EAF <- ifelse(BrS$EAF > 0.5, 1 - BrS$EAF, BrS$EAF)
setnames(BrS, old = 'EAF', new = 'MAF')

#save clean file in txt
write.table(BrS, "../clean/Barc_2022_BrS_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
