require(GenomicSEM)
require(data.table)
require(tidyverse)
library(tidyr)



data <- data.frame(
  samp.prevelance = c(0.08167744, 0.219, 0.123, 0.056, 0.316, 0.07981),
  pop.prevelance = c(0.017, 0.0001, 0.036, 0.092,0.22, 0.05),
  traitnames = c("HF", "BrS", "CAD", "STR", "HYP", "ANX"),
  trait_path = c("munge/HF.sumstats.gz", "munge/BrS.sumstats.gz", "munge/CAD.sumstats.gz", "munge/STR.sumstats.gz", "munge/HYP_2019.sumstats", "munge/ANX.sumstats.gz"),
  files = c("clean/Henry_2025_HF_EUR_clean.txt",
            "clean/Barc_2022_BrS_clean.txt",
            "clean/aragam_2022_CAD_clean.txt",
            "clean/mishra_2022_stroke_clean.txt",
            "clean/Zhu_2019_HYP_clean.txt",
            "clean/Friligkou_2024_ANX_EUR_clean.txt")
)

#ldsc
ld <- "eur_w_ld_chr/eur_w_ld_chr"
wld <- "eur_w_ld_chr/eur_w_ld_chr"

#LDSCoutput <- ldsc(data$trait_path, sample.prev = data$samp.prevelance, population.prev = data$pop.prevelance, ld, wld, data$traitnames)
#save(LDSCoutput, file = "LDSCoutput.Rdata")
load("LDSCoutput.Rdata")


#model all
model_all <- 'F1 =~ NA*CAD + STR + HF + BrS
F2 =~ A*F1 + A*ANX'

run_model_all <- usermodel(LDSCoutput, estimation = "DWLS", model = model_all, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

# other model seperate
model_seperate <- 'F1 =~ NA*CAD + STR + HF
F2 =~ A*BrS + A*HF 
F2 ~~ F1
F3 =~ B*F1 + B*ANX
F4 =~ C*F2 + C*ANX'

run_model_seperate <- usermodel(LDSCoutput, estimation = "DWLS", model = model_seperate, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

# model for userGWAS, std.lv = FALSE

model_all <- 'F1 =~ NA*CAD + STR + HF + BrS + HYP
          F2 =~ 0.27054220*F1 + 0.27054220*ANX
F1 ~~ 1*F1
F2 ~~ 1*F2'

model_seperate <- 'F1 =~ NA*CAD + STR + HF
F2 =~ A*BrS + A*HF 
F2 ~~ F1
F3 =~ 0.86047998*F1 + 0.86047998*ANX
F4 =~ 0.93108622*F2 + 0.93108622*ANX
F3 ~~ 1*F3
F4 ~~ 1*F4'





