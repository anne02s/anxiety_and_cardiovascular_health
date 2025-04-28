require(data.table)
library(tidyr)
require(tidyverse)

#load file
PR <- fread("../data/weng_2020_PR_unclean.txt")

#clean file:

#delete columns that are not necessary
print(paste("PR table delete columns:", nrow(PR)))
PR <- PR %>% select(-chr, -pos)
print(paste("PR table delete columns(no difference!):", nrow(PR)))

#rename columns
print(paste("PR table rename:", nrow(PR)))
colnames(PR) <- c("SNP", "A1", "A2", "EAF", "beta", "SE", "p_value") 
print(paste("PR table rename(no difference!):", nrow(PR)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(PR$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  PR$A1_multiple <- nchar(PR$A1) > 1
  # delete those SNPs
  PR <- PR[PR$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("PR table delete A1_multiple:", nrow(PR)))
  PR <- PR %>% select(-A1_multiple)
  print(paste("PR table delete A1_multiple(no difference!):", nrow(PR)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(PR$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  PR$A2_multiple <- nchar(PR$A2) > 1
  # delete those SNPs
  PR <- PR[PR$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("PR table delete A2_multiple:", nrow(PR)))
  PR <- PR %>% select(-A2_multiple)
  print(paste("PR table delete A2_multiple(no difference!):", nrow(PR)))
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- PR$SNP == "" | !grepl("^r", PR$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the once that are empty
  PR <- PR[!not_r_or_empty, ]
}

# Delete SNPs that are duplicates
PR <- PR[!duplicated(PR$SNP), ]

#make MAF
PR$EAF <- ifelse(PR$EAF > 0.5, 1 - PR$EAF, PR$EAF)
setnames(PR, old = 'EAF', new = 'MAF')

#calculate effective sample size implied by GWAS summary statistics 
PR$Neff<-4/((2*PR$MAF*(1-PR$MAF))*PR$SE^2)

#save clean file in txt
write.table(PR, "../clean/weng_2020_PR_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

