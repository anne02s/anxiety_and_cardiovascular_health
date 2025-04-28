require(data.table)
library(tidyr)
require(tidyverse)

#load file
SBP <- fread("../data/QC_SBP-meta-analysis_ICBP2024.tsv")

#clean file:

#rename columns
print(paste("SBP table rename:", nrow(SBP)))
colnames(SBP) <- c("chromosome", "bpl", "A1", "A2", "beta", "SE", "EAF", "p_value", "SNP", "N") 
print(paste("SBP table rename(no difference!):", nrow(SBP)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(SBP$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  SBP$A1_multiple <- nchar(SBP$A1) > 1
  # delete those SNPs
  SBP <- SBP[SBP$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("SBP table delete A1_multiple:", nrow(SBP)))
  SBP <- SBP %>% select(-A1_multiple)
  print(paste("SBP table delete A1_multiple(no difference!):", nrow(SBP)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(SBP$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  SBP$A2_multiple <- nchar(SBP$A2) > 1
  # delete those SNPs
  SBP <- SBP[SBP$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("SBP table delete A2_multiple:", nrow(SBP)))
  SBP <- SBP %>% select(-A2_multiple)
  print(paste("SBP table delete A2_multiple(no difference!):", nrow(SBP)))
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- SBP$SNP == "" | !grepl("^r", SBP$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the once that are empty
  print(paste("SBP table bevore delete empty:", nrow(SBP)))
  SBP <- SBP[!not_r_or_empty, ]
  print(paste("SBP table after delete empty:", nrow(SBP)))
}

# Delete SNPs that are duplicates
print(paste("SBP table bevore delete duplicates:", nrow(SBP)))
SBP <- SBP[!duplicated(SBP$SNP), ]
print(paste("SBP table after delete duplicated:", nrow(SBP)))

#make MAF
SBP$EAF <- ifelse(SBP$EAF > 0.5, 1 - SBP$EAF, SBP$EAF)
setnames(SBP, old = 'EAF', new = 'MAF')


#delete columns that are not necessary
print(paste("SBP table delete columns:", nrow(SBP)))
SBP <- SBP %>% select(-chromosome, -bpl)
print(paste("SBP table delete columns(no difference!):", nrow(SBP)))

#reorder columns
print(paste("SBP table reorder:", nrow(SBP)))
SBP <- SBP[, c('SNP', 'A1', 'A2', 'beta', 'SE', 'p_value', 'MAF', 'N')]
print(paste("SBP table reorder:(no difference!)", nrow(SBP)))

#save clean file in txt
write.table(SBP, "../clean/Keaton_2024_SBP_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
