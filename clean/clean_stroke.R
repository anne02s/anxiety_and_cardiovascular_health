require(data.table)
library(tidyr)
require(tidyverse)


# load file
STR <- fread("../data/mishra_2022_stroke_unclean.tsv_2")
print(paste("STR table:", nrow(STR)))
#clean file:

# rename column
print(paste("STR table rename:", nrow(STR)))
setnames(STR, old = 'effect_allele_frequency', new = 'EAF')
setnames(STR, old = 'standard_error', new = 'SE')
setnames(STR, old = 'effect_allele', new = 'A1')
setnames(STR, old = 'other_allele', new = 'A2')
print(paste("STR table rename(no difference!):", nrow(STR)))

# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(STR$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  STR$A1_multiple <- nchar(STR$A1) > 1
  # delete those SNPs
  STR <- STR[STR$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("STR table delete A1_multiple:", nrow(STR)))
  STR <- STR %>% select(-A1_multiple)
  print(paste("STR table delete A1_multiple(no difference!):", nrow(STR)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(STR$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  STR$A2_multiple <- nchar(STR$A2) > 1
  # delete those SNPs
  STR <- STR[STR$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("STR table delete A2_multiple:", nrow(STR)))
  STR <- STR %>% select(-A2_multiple)
  print(paste("STR table delete A2_multiple(no difference!):", nrow(STR)))
}

#make MAF
STR$EAF <- ifelse(STR$EAF > 0.5, 1 - STR$EAF, STR$EAF)
setnames(STR, old = 'EAF', new = 'MAF')

#delete columns that are not necessary
print(paste("STR table delete columns:", nrow(STR)))
STR <- STR %>% select(-ci_upper, -ci_lower, -odds_ratio)
print(paste("STR table delete columns(no difference!):", nrow(STR)))

#save clean file in txt
write.table(STR, "../clean/mishra_2022_stroke_no_IDrs_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



