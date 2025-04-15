require(data.table)
library(tidyr)
require(tidyverse)


# Coronary Artery Disease

#load file
CAD <- fread("../data/aragam_2022_CAD_unclean")
print(paste("CAD table:", nrow(CAD)))

#clean file:

# rename column
print(paste("CAD table rename:", nrow(CAD)))
setnames(CAD, old = 'rsID', new = 'SNP')
print(paste("CAD table rename(no difference!):", nrow(CAD)))

#reorder columns
print(paste("CAD table reorder:", nrow(CAD)))
CAD <- CAD[, c('SNP', 'A1', 'A2', 'beta', 'SE', 'p_value', 'N', 'EAF', 'bpchr', 'chromosome', 'base_pair_location')]
print(paste("CAD table reorder:(no difference!)", nrow(CAD)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(CAD$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  CAD$A1_multiple <- nchar(CAD$A1) > 1
  # delete those SNPs
  CAD <- CAD[CAD$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("CAD table delete A1_multiple:", nrow(CAD)))
  CAD <- CAD %>% select(-A1_multiple)
  print(paste("CAD table delete A1_multiple(no difference!):", nrow(CAD)))
}
  
# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(CAD$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  CAD$A2_multiple <- nchar(CAD$A2) > 1
  # delete those SNPs
  CAD <- CAD[CAD$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("CAD table delete A2_multiple:", nrow(CAD)))
  CAD <- CAD %>% select(-A2_multiple)
  print(paste("CAD table delete A2_multiple(no difference!):", nrow(CAD)))
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- CAD$SNP == "" | !grepl("^r", CAD$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the once that are empty
  CAD <- CAD[!not_r_or_empty, ]
}

# Delete SNPs that are duplicates
print(paste("CAD bevore deleting duplicates:", nrow(CAD)))
CAD <- CAD[!duplicated(CAD$SNP), ]
print(paste("CAD after deleting duplicates:", nrow(CAD)))

#make MAF
CAD$EAF <- ifelse(CAD$EAF > 0.5, 1 - CAD$EAF, CAD$EAF)
setnames(CAD, old = 'EAF', new = 'MAF')

#delete columns that are not necessary
print(paste("CAD table delete columns:", nrow(CAD)))
CAD <- CAD %>% select(-bpchr, -chromosome, -base_pair_location)
print(paste("CAD table delete columns(no difference!):", nrow(CAD)))


#save clean file in txt
write.table(CAD, "../clean/aragam_2022_CAD_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
