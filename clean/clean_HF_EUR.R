require(data.table)
library(tidyr)
require(tidyverse)

# load file
HF <- fread("../data/FORMAT-METAL_Pheno1_EUR.tsv")
print(paste("HF table:", nrow(HF)))
#clean file:

# rename column
print(paste("HF table rename:", nrow(HF)))
setnames(HF, old = 'rsID', new = 'SNP')
setnames(HF, old = 'A1_beta', new = 'BETA')
setnames(HF, old = 'pval', new = 'p_value')
setnames(HF, old = 'N_case', new = 'N')
setnames(HF, old = 'A1_freq', new = 'EAF')
setnames(HF, old = 'se', new = 'SE')
print(paste("HF table rename(no difference!):", nrow(HF)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(HF$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  HF$A1_multiple <- nchar(HF$A1) > 1
  # delete those SNPs
  HF <- HF[HF$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("HF table delete A1_multiple:", nrow(HF)))
  HF <- HF %>% select(-A1_multiple)
  print(paste("HF table delete A1_multiple(no difference!):", nrow(HF)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(HF$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  HF$A2_multiple <- nchar(HF$A2) > 1
  # delete those SNPs
  HF <- HF[HF$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("HF table delete A2_multiple:", nrow(HF)))
  HF <- HF %>% select(-A2_multiple)
  print(paste("HF table delete A2_multiple(no difference!):", nrow(HF)))
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- HF$SNP == "" | !grepl("^r", HF$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the empty once
  HF <- HF[!not_r_or_empty, ]
}

# Delete SNPs that are duplicates
print(paste("HF bevore deleting duplicates:", nrow(HF)))
HF <- HF[!duplicated(HF$SNP), ]
print(paste("HF after deleting duplicates:", nrow(HF)))

#make MAF
HF$EAF <- ifelse(HF$EAF > 0.5, 1 - HF$EAF, HF$EAF)
setnames(HF, old = 'EAF', new = 'MAF')

#delete columns that are not necessary
print(paste("HF table delete columns:", nrow(HF)))
HF <- HF %>% select(-"#key", -logP, -N_total, -isq_het, -p_het, -chr, -pos_b37)
print(paste("HF table delete columns(no difference!):", nrow(HF)))

#reorder columns
print(paste("HF table reorder:", nrow(HF)))
HF <- HF[, c('SNP', 'A1', 'A2', 'BETA', 'SE', 'p_value', 'N', 'MAF')]
print(paste("HF table reorder:(no difference!)", nrow(HF)))


#save clean file in txt
write.table(HF, "../clean/Henry_2025_HF_EUR_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
