require(data.table)
library(tidyr)
require(tidyverse)

#load file
DBP <- fread("../data/QC_DBP-meta-analysis_ICBP2024.tsv")

#clean file:

#rename columns
print(paste("DBP table rename:", nrow(DBP)))
colnames(DBP) <- c("chromosome", "bpl", "A1", "A2", "beta", "SE", "EAF", "p_value", "SNP", "N") 
print(paste("DBP table rename(no difference!):", nrow(DBP)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(DBP$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  DBP$A1_multiple <- nchar(DBP$A1) > 1
  # delete those SNPs
  DBP <- DBP[DBP$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("DBP table delete A1_multiple:", nrow(DBP)))
  DBP <- DBP %>% select(-A1_multiple)
  print(paste("DBP table delete A1_multiple(no difference!):", nrow(DBP)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(DBP$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  DBP$A2_multiple <- nchar(DBP$A2) > 1
  # delete those SNPs
  DBP <- DBP[DBP$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("DBP table delete A2_multiple:", nrow(DBP)))
  DBP <- DBP %>% select(-A2_multiple)
  print(paste("DBP table delete A2_multiple(no difference!):", nrow(DBP)))
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- DBP$SNP == "" | !grepl("^r", DBP$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the once that are empty
  print(paste("DBP table bevore delete empty:", nrow(DBP)))
  DBP <- DBP[!not_r_or_empty, ]
  print(paste("DBP table after delete empty:", nrow(DBP)))
}

# Delete SNPs that are duplicates
print(paste("DBP table bevore delete duplicates:", nrow(DBP)))
DBP <- DBP[!duplicated(DBP$SNP), ]
print(paste("DBP table after delete duplicates:", nrow(DBP)))

#make MAF
DBP$EAF <- ifelse(DBP$EAF > 0.5, 1 - DBP$EAF, DBP$EAF)
setnames(DBP, old = 'EAF', new = 'MAF')


#delete columns that are not necessary
print(paste("DBP table delete columns:", nrow(DBP)))
DBP <- DBP %>% select(-chromosome, -bpl)
print(paste("DBP table delete columns(no difference!):", nrow(DBP)))

#reorder columns
print(paste("DBP table reorder:", nrow(DBP)))
DBP <- DBP[, c('SNP', 'A1', 'A2', 'beta', 'SE', 'p_value', 'MAF', 'N')]
print(paste("DBP table reorder:(no difference!)", nrow(DBP)))

#save clean file in txt
write.table(DBP, "../clean/Keaton_2024_DBP_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
