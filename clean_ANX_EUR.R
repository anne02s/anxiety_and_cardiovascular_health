require(data.table)
library(tidyr)
require(tidyverse)



# during the cleaning we are goining to do these things:
# *: check of row count does not change during alterations where change is not needed
# *: delete unnecessary columns
# *: check is all SNP are single nucleotide #(any(nchar(file$A1) > 1))
# *: check if all SNP start with rs #(any(!grepl("^r", CAD$SNP, ignore.case = TRUE)))
# *: make new file with clean data

# load file
ANX_EUR <- fread("../data/ANX_EUR.txt")
print(paste("ANX_EUR table:", nrow(ANX_EUR)))
#clean file:

#reorder columns
print(paste("ANX_EUR table reorder:", nrow(ANX_EUR)))
ANX_EUR <- ANX_EUR[, c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'Z', 'MAF', 'CHR', 'BP')]
print(paste("ANX_EUR table reorder:(no difference!)", nrow(ANX_EUR)))


# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(ANX_EUR$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  ANX_EUR$A1_multiple <- nchar(ANX_EUR$A1) > 1
  # delete those SNPs
  ANX_EUR <- ANX_EUR[ANX_EUR$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("ANX_EUR table delete A1_multiple:", nrow(ANX_EUR)))
  ANX_EUR <- ANX_EUR %>% select(-A1_multiple)
  print(paste("ANX_EUR table delete A1_multiple(no difference!):", nrow(ANX_EUR)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(ANX_EUR$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  ANX_EUR$A2_multiple <- nchar(ANX_EUR$A2) > 1
  # delete those SNPs
  ANX_EUR <- ANX_EUR[ANX_EUR$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("ANX_EUR table delete A2_multiple:", nrow(ANX_EUR)))
  ANX_EUR <- ANX_EUR %>% select(-A2_multiple)
  print(paste("ANX_EUR table delete A2_multiple(no difference!):", nrow(ANX_EUR)))
}

#calculate effective sample size implied by GWAS summary statistics 
ANX_EUR$Neff <- 4/((2*ANX_EUR$MAF*(1-ANX_EUR$MAF))*ANX_EUR$SE^2)

# Delete SNPs that are duplicates
ANX_EUR <- ANX_EUR[!duplicated(ANX_EUR$SNP), ]

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- ANX_EUR$SNP == "" | !grepl("^r", ANX_EUR$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
}

#delete columns that are not necessary
print(paste("ANX_EUR table delete columns:", nrow(ANX_EUR)))
ANX_EUR <- ANX_EUR %>% select(-CHR, -BP)
print(paste("ANX_EUR table delete columns(no difference!):", nrow(ANX_EUR)))

#save clean file in txt
write.table(ANX_EUR, "../clean/Friligkou_2024_ANX_EUR_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
