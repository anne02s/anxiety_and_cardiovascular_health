library(dplyr)
require(tidyverse)
library(tidyr)
require(data.table)

# load the data
gwas_data <- fread("../clean/mishra_2022_stroke_no_IDrs_clean.txt")
reference_panel <- fread("../reference.1000G.maf.0.005.txt")

# rename the colums of the reference panel so it is easy to merge
print(paste("reference_panel rename:", nrow(reference_panel)))
colnames(reference_panel) <- c("SNP", "chromosome", "base_pair_location", "M_A_F", "A_1", "A_2")
print(paste("reference_panel rename(no difference!):", nrow(reference_panel)))

# merge the data
print(paste("gwas_data merge:", nrow(gwas_data)))
merged_data <- left_join(gwas_data, reference_panel, by = c("chromosome", "base_pair_location"))
print(paste("merged_data:", nrow(merged_data)))


# Check if there are any non-NA mismatches between A1 and A_1
if (any(merged_data$A1 != merged_data$A_1, na.rm = TRUE)) {
  # Mark rows where A1 and A_1 are not identical AND A_1 is not NA
  merged_data$A1_not_identical <- !is.na(merged_data$A_1) & merged_data$A1 != merged_data$A_1

  # Replace the values in A1 with the values from A_1 where they are not identical
  merged_data$A1[merged_data$A1_not_identical] <- merged_data$A_1[merged_data$A1_not_identical]
}


# Check if there are any non-NA mismatches between A2 and A_2
if (any(merged_data$A2 != merged_data$A_2, na.rm = TRUE)) {
  # Mark rows where A2 and A_2 are not identical AND A_2 is not NA
  merged_data$A2_not_identical <- !is.na(merged_data$A_2) & merged_data$A2 != merged_data$A_2

  # Replace the values in A2 with the values from A_2 where they are not identical
  merged_data$A2[merged_data$A2_not_identical] <- merged_data$A_2[merged_data$A2_not_identical]
}

# delete beta if it does not exsist
print(paste("table bevore deleting NA in beta:", nrow(merged_data)))
merged_data <- merged_data[!is.na(merged_data$beta), ]
print(paste("table after deleting NA in beta:", nrow(merged_data)))

# flip beta and MAF if A1 and A2 are changed
if (any(merged_data$A1_not_identical == TRUE & merged_data$A2_not_identical == TRUE)) {
	count <- sum(merged_data$A1_not_identical == TRUE & merged_data$A2_not_identical == TRUE)
        print(paste("count change beta and MAF:", count))
	# flip beta and MAF
        idx <- which(merged_data$A1_not_identical & merged_data$A2_not_identical)
	merged_data$beta[idx] <- -merged_data$beta[idx]
        merged_data$MAF[idx] <- 1 - merged_data$MAF[idx]
}

# Check if there are SNPs that are empty or do not start with 'r'
not_r_or_empty <- merged_data$SNP == "" | !grepl("^r", merged_data$SNP, ignore.case = TRUE)

if (any(not_r_or_empty)) {
  SNPs <- sum(not_r_or_empty)
  print(SNPs)
  # delete the once that are empty
  merged_data <- merged_data[!not_r_or_empty, ]
}

# Delete SNPs that are duplicates
print(paste("merged_data bevore deleting duplicates:", nrow(merged_data)))
merged_data <- merged_data[!duplicated(merged_data$SNP), ]
print(paste("merged_data after deleting duplicates:", nrow(merged_data)))

#calculate effective sample size implied by GWAS summary statistics 
STR$N<-4/((2*STR$MAF*(1-STR$MAF))*STR$SE^2)

#delete the columns 
print(paste("merged_data delete column:", nrow(merged_data)))
merged_data <- merged_data %>% select(-A1_not_identical, -A2_not_identical, -chromosome, -base_pair_location, -M_A_F)
print(paste("merged_data delete column(no difference!):", nrow(merged_data)))

#reorder columns
print(paste("merged_data table reorder:", nrow(merged_data)))
merged_data <- merged_data[, c('SNP', 'A1', 'A2', 'beta', 'SE', 'p_value', 'MAF', 'N')]
print(paste("merged_data table reorder:(no difference!)", nrow(merged_data)))

#write table
write.table(merged_data, "../clean/mishra_2022_stroke_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
