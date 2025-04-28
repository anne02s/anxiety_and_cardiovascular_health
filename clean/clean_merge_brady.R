library(dplyr)
require(tidyverse)
library(tidyr)
require(data.table)

# load the data
#bra_dist_inc <- fread("../data/bradyarrhythmia/DISTAL_INC_filtered.txt")
#bra_dist_rest <- fread("../data/bradyarrhythmia/DISTAL_REST_filtered.txt")
#bra_pacer <- fread("../data/bradyarrhythmia/PACER_filtered.txt")
#bra_snd_inc <- fread("../data/bradyarrhythmia/SND_INC_filtered.txt")
bra_snd_rest <- fread("../data/bradyarrhythmia/SND_REST_filtered.txt")
reference_panel <- fread("../reference.1000G.maf.0.005.txt")


# seperate the fist colum 
print(paste("seprate colum:", nrow(bra_snd_rest)))
gwas_data <- bra_snd_rest %>% tidyr::separate(MarkerName, into = c("chromosome", "base_pair_location"), sep = ":")
print(paste("seprate colum:(no difference!)", nrow(gwas_data)))

# rename the colums 
print(paste("table rename:", nrow(gwas_data)))
setnames(gwas_data, old = c("Effect", "StdErr", "Allele1", "Allele2"),
                new = c("beta", "SE", "A1", "A2"))
print(paste("table rename:(no difference!)", nrow(gwas_data)))

# rename the colums of the reference panel so it is easy to merge
print(paste("reference_panel rename:", nrow(reference_panel)))
colnames(reference_panel) <- c("SNP", "chromosome", "base_pair_location", "MAF", "A_1", "A_2")
print(paste("reference_panel rename(no difference!):", nrow(reference_panel)))
    
#write table
write.table(gwas_data, "../clean/weng_2025_bra_snd_rest.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
gwas_data_ <- fread("../clean/weng_2025_bra_snd_rest.txt")

# in brady they added zero's to make sure all are the same legth. this gives an error.
reference_panel$base_pair_location <- sprintf("%09d", as.numeric(reference_panel$base_pair_location))
gwas_data_$base_pair_location <- sprintf("%09d", as.numeric(gwas_data_$base_pair_location))

# Just in case fread messed up:
gwas_data_$chromosome <- as.character(gwas_data_$chromosome)
reference_panel$chromosome <- as.character(reference_panel$chromosome)

# merge the data
print(paste("gwas_data merge:", nrow(gwas_data_)))
merged_data <- left_join(gwas_data_, reference_panel, by = c("chromosome", "base_pair_location"))
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

# delete if beta does not exsist
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
merged_data$N<-4/((2*merged_data$MAF*(1-merged_data$MAF))*merged_data$SE^2)

#delete the columns 
print(paste("merged_data delete column:", nrow(merged_data)))
merged_data <- merged_data %>% select(-A1_not_identical, -A2_not_identical, -chromosome, -base_pair_location)
print(paste("merged_data delete column(no difference!):", nrow(merged_data)))

#reorder columns
print(paste("merged_data table reorder:", nrow(merged_data)))
merged_data <- merged_data[, c('SNP', 'A1', 'A2', 'beta', 'SE', 'P-value', 'MAF', 'N')]
print(paste("merged_data table reorder:(no difference!)", nrow(merged_data)))

#write table
write.table(merged_data, "../clean/weng_2025_bra_snd_rest_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


