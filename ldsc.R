require(GenomicSEM)
require(data.table)
require(tidyverse)
library(tidyr)

HF_prev <- 0.08167744
AF_prev <- 166322/1480272
brady_dist_inc_prev <- 28359/812226
brady_dist_rest_prev <- 7834/1112609
brady_snd_inc_prev <- 8891/1117885
brady_snd_rest_prev <- 4622/1108168
brady_pace_prev <- 26885/1142609
BrS_prev <- 2820/12821
ANX_prev <- 0.07981
CAD_prev <- 120788/980319
STR_prev <- 73652/1308460
HYP_prev <- 0.316
QRS_prev <- NA
QT_prev <- NA


data <- data.frame(
  samp.prevelance = c(HF_prev, AF_prev, brady_dist_inc_prev, brady_dist_rest_prev, brady_pace_prev, brady_snd_inc_prev, brady_snd_rest_prev, BrS_prev, CAD_prev, STR_prev, ANX_prev, HYP_prev, QRS_prev, QT_prev),
  pop.prevelance = c(0.017, 0.02, 0.038, 0.038, 0.038, 0.038, 0.038, 0.0001, 0.036, 0.092, 0.05, 0.22, NA, NA),
  traitnames <- c("HF", "AF", "BRA_dist_inc", "BRA_dist_rest", "BRA_pacer", "BRA_snd_inc", "BRA_snd_rest", "BrS", "CAD", "STR", "ANX","HYP", "QRS". "QT"),
  trait_path <- c("munge/HF.sumstats.gz", "munge/AF.sumstats.gz", "munge/brady_dist_inc.sumstats.gz", "munge/brady_dist_rest.sumstats.gz", "munge/brady_pacer.sumstats.gz", "munge/brady_snd_inc.sumstats.gz", "munge/brady_snd_rest.sumstats.gz","munge/BrS.sumstats.gz", "munge/CAD.sumstats.gz", "munge/STR.sumstats.gz", "munge/ANX.sumstats.gz", "munge/HYP.sumstats.gz", "munge/QRS.sumstats.gz", "munge/QT.sumstats.gz")
  
)



ld <- "eur_w_ld_chr/eur_w_ld_chr"
wld <- "eur_w_ld_chr/eur_w_ld_chr"

LDSCoutput <- ldsc(data$trait_path, sample.prev = data$samp.prevelance, population.prev = data$pop.prevelance, ld, wld, data$traitnames)

save(LDSCoutput, file = "LDSC_output.Rdata")

#heatmap
library(ggplot2)
library(reshape2)
mat <- cov2cor(LDSCoutput$S)
data_melted <- as.data.frame(as.table(mat))
data_melted$Var1 <- data$traitnames
data_melted$Var1 <- factor(data_melted$Var1, levels = unique(data_melted$Var1))
data_melted$Var2 <- factor(data_melted$Var2, levels = unique(data_melted$Var2))
# Create a heatmap using ggplot2
ggplot(data_melted, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Freq)), color = "black", size = 6, vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0.10,
    limits = c(min(data_melted$Freq), max(data_melted$Freq))
  ) +
  theme_minimal() +
  labs(title = "cov2cor", x = "X-axis", y = "Y-axis")


#calculating standard error: 
n<-nrow(LDSCoutput$S)
SE<-matrix(0, n, n)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))

# Calculate z-scores
z <- LDSCoutput$S / SE

#heatmap
z[z == Inf] <- NA
Z <- melt(z)
Z$Var1 <- data$traitnames
Z$Var1 <- factor(Z$Var1, levels = unique(Z$Var1))
Z$Var2 <- factor(Z$Var2, levels = unique(Z$Var2))
# Clean up your data first
Z_clean <- Z
Z_clean$value[is.infinite(Z_clean$value)] <- NA  # Turn -Inf to NA
Z_clean$label <- ifelse(is.na(Z_clean$value), "", sprintf("%.2f", Z_clean$value))

ggplot(Z_clean, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black", size = 6, vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 4.00,
    na.value = "white"  # This ensures NA or -Inf are white tiles
  ) +
  scale_y_discrete(position = "right") +  
  theme_minimal() +
  labs(title = "Z", x = "X-axis", y = "Y-axis")

# Calculate two-tailed p-values
p <- 2 * (1 - pnorm(abs(z)))
P <- melt(p)
P$Var1 <- data$traitnames
P$Var1 <- factor(P$Var1, levels = unique(P$Var1))
P$Var2 <- factor(P$Var2, levels = unique(P$Var2))
# Create a heatmap using ggplot2
# Clean up your data first
P_clean <- P
P_clean$value[is.infinite(P_clean$value)] <- NA  # Turn -Inf to NA
P_clean$label <- ifelse(is.na(P_clean$value), "", sprintf("%.2f", P_clean$value))

ggplot(P_clean, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black", size = 6, vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(
    low = "red",
    mid = "white",
    high = "blue",
    midpoint = 0.05,
    na.value = "white"  # This ensures NA or -Inf are white tiles
  ) +
  scale_y_discrete(position = "right") +  
  theme_minimal() +
  labs(title = "P", x = "X-axis", y = "Y-axis")
