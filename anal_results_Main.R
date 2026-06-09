rm(list=ls())

library(dplyr)
library(lmtp)
library(mice)
library(glmnet)
library(purrr)

# part0 - load data ----
rdata_files <- list.files(
  pattern = "\\.RData$",
  full.names = TRUE
)

#for lmtp
out <- vector("list", length(rdata_files))
summary_df <- NULL

#for effect modificatop
pooled_results <- list()

for (f in rdata_files) {
  
  load(f)
  
#load("Tmin_CF3_trim0.9_trt_suffix_count_above_95.RData")
# results_lmtp[[1]]$contrast$vals$theta
# results_lmtp[[1]]$contrast$vals$std.error

# results_lmtp[[1]]$contrast$eifs %>% head
# (results_lmtp[[1]]$factual$eif - results_lmtp[[1]]$cfactual$eif) %>% head


# part I - pooled effect estimates ----
# Number of imputations
m <- 10
Exp <- "Tmin"#Tmean, Tmax, Tmin
#data_path <- "C://Users/Guillaume/Desktop/EnvEpi/BAS-CDI-SDQ/BAS/"
data_path <- "/bettik/barbalag/BAS/"


# Extract theta and SE
theta <- sapply(1:m, function(i) {
  results_lmtp[[i]]$contrast$vals$theta
})

se <- sapply(1:m, function(i) {
  results_lmtp[[i]]$contrast$vals$std.error
})

# Rubin's rules
theta_pooled <- mean(theta)

U_bar <- mean(se^2)
B <- var(theta)

T_var <- U_bar + (1 + 1/m) * B
se_pooled <- sqrt(T_var)

# Results
list(
  theta_pooled = theta_pooled,
  se_pooled    = se_pooled,
  m            = m,
  U_bar        = U_bar,
  B            = B,
  total_var    = T_var
)

lambda <- ((1 + 1/m) * B) / T_var
df <- (m - 1) / lambda^2
alpha <- 0.05
crit <- qt(1 - alpha/2, df)
CI <- c(
  lower = theta_pooled - crit * se_pooled,
  upper = theta_pooled + crit * se_pooled
)

# Test statistic
t_stat <- theta_pooled / se_pooled
# Two-sided p-value
p_value <- 2 * (1 - pt(abs(t_stat), df))

out[[f]] <- tibble(
  file = basename(f),
  theta_pooled = theta_pooled,
  se_pooled    = se_pooled,
  CI_low       = CI[1],
  CI_high      = CI[2],
  p_value      = p_value
)

# partII - effect modification see anal_results_effect_mdif.R ----

summary_df <- bind_rows(out)
summary_df
save(summary_df, file="BAS_lmtp_Main.RData")