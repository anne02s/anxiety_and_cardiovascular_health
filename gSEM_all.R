require(GenomicSEM)
require(data.table)
require(tidyverse)
library(tidyr)



data <- data.frame(
  samp.prevelance = c(0.08167744, 2820/12821, 120788/980319, 73652/1308460, 0.07981),
  pop.prevelance = c(0.0172, 0.05, 0.063, 0.019, 0.05),
  traitnames = c("HF", "BrS", "CAD", "STR", "ANX"),
  trait_path = c("munge/HF.sumstats.gz", "munge/BrS.sumstats.gz", "munge/CAD.sumstats.gz", "munge/STR.sumstats.gz", "munge/ANX.sumstats.gz"),
  files = c("clean/Henry_2025_HF_EUR_clean.txt",
            "clean/Barc_2022_BrS_clean.txt",
            "clean/aragam_2022_CAD_clean.txt",
            "clean/mishra_2022_stroke_clean.txt",
             "clean/Friligkou_2024_ANX_EUR_clean.txt")
)


#ldsc
ld <- "eur_w_ld_chr/eur_w_ld_chr"
wld <- "eur_w_ld_chr/eur_w_ld_chr"

#LDSCoutput <- ldsc(data$trait_path, sample.prev = data$samp.prevelance, population.prev = data$pop.prevelance, ld, wld, data$traitnames)
#save(LDSCoutput, file = "LDSCoutput.Rdata")
load("LDSCoutput.Rdata")



#model all
model <- 'F1 =~ NA*CAD + STR + HF + BrS
F2 =~ A*F1 + A*ANX
F2 =~ 1*F2'

run_model <- usermodel(LDSCoutput, estimation = "DWLS", model = model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
run_model



# other model seperate
model <- 'F1 =~ NA*CAD + STR + HF
F2 =~ A*BrS + A*HF 
F2 ~~ F1
F3 =~ B*F1 + B*ANX
F4 =~ C*F2 + C*ANX
F3 ~~ 1*F3
F4 ~~ 1*F4'

model <- 'F1 =~ NA*CAD + STR + HF
F2 =~ A*BrS + A*HF 
F2 ~~ F1
F3 =~ B*F1 + B*ANX
F4 =~ C*F2 + C*ANX
F5 =~ D*F4 + D*F3
F5 ~~ 1*F5'

           
ref= "reference.1000G.maf.0.005.txt"
se.logit=c(T,T,T,T,T)
Hail=c(F,F,F,F,F)
info.filter=.6
maf.filter=0.01
OLS=NULL
N=NULL
betas=NULL


PSYCH_sumstats <-sumstats(files=data$files,ref=ref,trait.names=data$traitnames,se.logit=se.logit,OLS=NULL,linprob=Hail,N=N,betas=NULL, info.filter=.6, maf.filter=0.01,keep.indel=FALSE,parallel=TRUE,cores=NULL)
save(PSYCH_sumstats, file = "PSYCH_sumstats.Rdata")
load("PSYCH_sumstats.Rdata")
#note that the code written above specifies all arguments for completeness, but as many of these arguments
#are set to the package default it could also be written as below and produce the same result:
#PSYCH_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = PSYCH_sumstats, estimation = "DWLS", parallel=TRUE)
PSYCH_factor <- userGWAS(covstruc = LDSCoutput,
                         SNPs = PSYCH_sumstats,
                         estimation = "DWLS",
                         model = model,
                         printwarn = TRUE,
                         cores = 1,
                         toler = FALSE,
                         SNPSE = FALSE,
                         parallel = FALSE,
                         fix_measurement=T)
PSYCH_facto




