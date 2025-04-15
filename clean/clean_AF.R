require(data.table)
library(tidyr)
require(tidyverse)

#load file
AF <- fread("../data/Miyazawa_2023_AF_unclean.tsv")

#clean file:

#rename columns
print(paste("AF table rename:", nrow(AF)))
colnames(AF) <- c("SNP", "chromosome", "bpl", "A1", "A2", "EAF", "beta", "SE", "p_value") 
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
}

# Delete SNPs that are duplicates
AF <- AF[!duplicated(AF$SNP), ]

#make MAF
AF$EAF <- ifelse(AF$EAF > 0.5, 1 - AF$EAF, AF$EAF)
setnames(AF, old = 'EAF', new = 'MAF')

#calculate effective sample size implied by GWAS summary statistics 
AF$Neff<-4/((2*AF$MAF*(1-AF$MAF))*AF$SE^2)

#delete columns that are not necessary
print(paste("AF table delete columns:", nrow(AF)))
AF <- AF %>% select(-chromosome, -bpl)
print(paste("AF table delete columns(no difference!):", nrow(AF)))

#save clean file in txt
write.table(AF, "../clean/Miyazawa_2023_AF_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
