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

rdata_files <- rdata_files[
  grepl("Tmin", rdata_files) &
  grepl("consec_count_above_95", rdata_files)
]

Dx_all <- list()
Dx_summary <- list()
totalN <- 8965

for (f in rdata_files) {
  
  print(f)
  load(f)

   Dx_all[[f]] <- map_dfr(seq_along(results_lmtp), function(i) {
    
    w <- results_lmtp[[i]]$cfactual$density_ratios %>% as.data.frame

    tibble(
      period = colnames(w),
      mean = apply(w, 2, function(x) {
        x <- (x[!is.na(x)]/mean(x))
        mean(x)
      }),
      median = apply(w, 2, function(x) {
        x <- (x[!is.na(x)]/mean(x))
        median(x)
      }),
      max = apply(w, 2, function(x) {
        x <- x[!is.na(x)]/mean(x)
        max(x)        
        }),
      min = apply(w, 2, function(x) {
        x <- x[!is.na(x)]/mean(x)
        min(x)
      }),
      
      MWP = apply(w, 2, function(x) {
        x <- x[!is.na(x)]/mean(x)
        max(x) / sum(x)
      }),
      CV = apply(w, 2, function(x) {
        x <- x[!is.na(x)]/mean(x)
        sd(x) / mean(x)
      }),
      ESS = apply(w, 2, function(x) {
        x <- x[!is.na(x)]/mean(x)
        (sum(x)^2) / sum(x^2)/totalN
      }),
      
      imputation = i
    )
  })
  
  Dx_summary[[f]] <- Dx_all[[f]] %>%
    group_by(period) %>%
    summarise(
      mean = mean(mean),
      median = mean(median),
      min = mean(min),
      max = mean(max),
      
      
      mean_ESS = mean(ESS),
      min_ESS  = min(ESS),
      max_ESS  = max(ESS),
      mean_CV = mean(CV),
      min_CV  = min(CV),
      max_CV  = max(CV),
      mean_MWP = mean(MWP),
      min_MWP  = min(MWP),
      max_MWP  = max(MWP),
    )
  print(Dx_summary[[f]])
}
save(Dx_all,Dx_summary,file="BAS_DR_Dx_all.RData")