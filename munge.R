require(GenomicSEM)
require(data.table)
require(tidyverse)
library(tidyr)


data <- data.frame(
  file = c("../clean/mishra_2022_stroke_clean.txt",
           "../clean/Friligkou_2024_ANX_EUR_clean.txt",
           "../clean/aragam_2022_CAD_clean.txt",
           "../clean/Henry_2025_HF_EUR_clean.txt"),
  traitnames = c("STR", "ANX", "CAD", "HF"),
  N = c(NA, NA, NA, NA),
  cases = c(73652, 87517, 120788, 153174)
)


munge(file = data$file , hm3 = "w_hm3.snplist", trait.names = data$traitnames, N = data$N, info.filter =  0.9, maf.filter =  0.01)