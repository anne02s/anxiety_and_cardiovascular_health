require(data.table)
library(tidyr)
require(tidyverse)

#load file
AF <- fread("../data/Roselli_2025_AF_unclean.txt")

#clean file:

# make p
AF$P <- exp(AF$`log(P)`)

#delete columns that are not necessary
print(paste("AF table delete columns:", nrow(AF)))
AF <- AF %>% select(-MarkerName, -chr, -position_b38, -`log(P)`, -n_events, -mean_impQual)
print(paste("AF table delete columns(no difference!):", nrow(AF)))


#rename columns
print(paste("AF table rename:", nrow(AF)))
colnames(AF) <- c("SNP", "A1", "A2", "MAF", "beta", "SE", "N", "p_value") 
print(paste("AF table rename(no difference!):", nrow(AF)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(AF$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  AF$A1_multiple <- nchar(AF$A1) > 1
  # delete those SNPs
  AF <- AF[AF$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("AF table delete A1_multiple:", nrow(AF)))
  AF <- AF %>% select(-A1_multiple)
  print(paste("AF table delete A1_multiple(no difference!):", nrow(AF)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(AF$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  AF$A2_multiple <- nchar(AF$A2) > 1
  # delete those SNPs
  AF <- AF[AF$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("AF table delete A2_multiple:", nrow(AF)))
  AF <- AF %>% select(-A2_multiple)
  print(paste("AF table delete A2_multiple(no difference!):", nrow(AF)))
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- AF$SNP == "" | !grepl("^r", AF$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the once that are empty
  AF <- AF[!not_r_or_empty, ]
}

# Delete SNPs that are duplicates
AF <- AF[!duplicated(AF$SNP), ]
print(paste("AF table after delete duplicated:", nrow(AF)))

#save clean file in txt
write.table(AF, "../clean/Roselli_2025_AF_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
