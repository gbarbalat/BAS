#changes
#/*
#*/

rm(list=ls())
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mice)
library(data.table)
library(stringr)
library(purrr)
library(ggplot2)

# header ----
file_ext <- ".sas7bdat"
Exp <- "Tmin"

load(paste0(Exp, "_exp.RData"))
Exp_preN_list <- get(paste0(Exp, "_exp_preN_list"))
Exp_postN_list <- get(paste0(Exp, "_exp_postN_list"))

# preproc step re- exposure ----


## split postN into 3 periods ----
# Function to extract a period
extract_period <- function(x, start, end) {
  n <- length(x)
  if (is.null(x) | start>n) return(NULL)  # no data in this range
  x[start:min(end, n)]
}

# Create new lists

#/*
Exp_preN_list1 <- lapply(Exp_preN_list, extract_period, start = 1, end = 13*7)
Exp_preN_list2 <- lapply(Exp_preN_list, extract_period, start = 13*7+1, end = 27*7)
Exp_preN_list3 <- lapply(Exp_preN_list, extract_period, start = 27*7+1, end = Inf)

names(Exp_preN_list1) <- names(Exp_preN_list)
names(Exp_preN_list2) <- names(Exp_preN_list)
names(Exp_preN_list3) <- names(Exp_preN_list)

length_Exp_preN3 <- sapply(Exp_preN_list, length) - (27*7)
#*/


Exp_postN_list1 <- lapply(Exp_postN_list, extract_period, start = 1, end = 365)
Exp_postN_list2 <- lapply(Exp_postN_list, extract_period, start = 366, end = 365+366)
Exp_postN_list3 <- lapply(Exp_postN_list, extract_period, start = 365+366+1, end = Inf)

# Preserve names
names(Exp_postN_list1) <- names(Exp_postN_list)
names(Exp_postN_list2) <- names(Exp_postN_list)
names(Exp_postN_list3) <- names(Exp_postN_list)

## Length of all elements in Exp ----
length_Exp_postN3 <- sapply(Exp_postN_list, length) - (365+366)


## save ----
#/*
save(length_Exp_preN3, length_Exp_postN3, file=paste0(Exp,"_length_Exp_3monthsPregnancy.RData"))
#*/

## Ndays above thresh ----
# Combine all values across all participants to compute global percentiles
all_values <- unlist(c(Exp_preN_list,Exp_postN_list))

p90 <- quantile(all_values, 0.9, na.rm = TRUE)
p95 <- quantile(all_values, 0.95, na.rm = TRUE)
p10 <- quantile(all_values, 0.1, na.rm = TRUE)
p05 <- quantile(all_values, 0.05, na.rm = TRUE)

above_thresh <- function(Exp) {
  
  # Function to count temperature above threshold
  count_above <- function(x, threshold) {
    if (is.null(x) | length(x) < 10) return(NA)  # no data in this range
    
    sum(x > threshold, na.rm = TRUE)
  }
  
  # Count consecutive days above threshold
  count_above_consecutive <- function(x, threshold) {
    if (is.null(x) | length(x) < 10) return(NA)  # no data in this range
    
    # Logical vector: TRUE if above threshold
    above <- x > threshold
    
    # Shift vector by 1 and check where both today and tomorrow are above
    sum(above[-length(above)] & above[-1], na.rm = TRUE)
  }
  
  # Apply to each element in the list
  count_90 <- sapply(Exp, count_above, threshold = p90)
  count_95 <- sapply(Exp, count_above, threshold = p95)
  consec_count_90 <- sapply(Exp, count_above_consecutive, threshold = p90)
  consec_count_95 <- sapply(Exp, count_above_consecutive, threshold = p95)
  
  # Combine in a data frame for clarity
  result <- data.frame(
    element = names(Exp),
    count_above_90 = count_90,
    count_above_95 = count_95,
    consec_count_above_90 = consec_count_90,
    consec_count_above_95 = consec_count_95
  )
  
  result
}
below_thresh <- function(Exp) {
  
  # Function to count temperature below threshold
  count_below <- function(x, threshold) {
    if (is.null(x) | length(x) < 10 ) return(NA)  # no data in this range
    
    sum(x < threshold, na.rm = TRUE)
  }
  
  # Count consecutive days below threshold
  count_below_consecutive <- function(x, threshold) {
    if (is.null(x) | length(x) < 10) return(NA)  # no data in this range
    
    # Logical vector: TRUE if below threshold
    below <- x < threshold
    
    # Shift vector by 1 and check where both today and tomorrow are below
    sum(below[-length(below)] & below[-1], na.rm = TRUE)
  }
  
  # Apply to each element in the list
  count_10 <- sapply(Exp, count_below, threshold = p10)
  count_05 <- sapply(Exp, count_below, threshold = p05)
  consec_count_10 <- sapply(Exp, count_below_consecutive, threshold = p10)
  consec_count_05 <- sapply(Exp, count_below_consecutive, threshold = p05)
  # count_90 <- sapply(Exp, count_above, threshold = p90)
  # count_95 <- sapply(Exp, count_above, threshold = p95)
  # consec_count_90 <- sapply(Exp, count_above_consecutive, threshold = p90)
  # consec_count_95 <- sapply(Exp, count_abpve_consecutive, threshold = p95)
  
  # Combine in a data frame for clarity
  result <- data.frame(
    element = names(Exp),
    count_below_10 = count_10,
    count_below_05 = count_05,
    consec_count_below_10 = consec_count_10,
    consec_count_below_05 = consec_count_05
    # count_above_90 = count_90,
    # count_above_95 = count_95,
    # consec_count_above_90 = consec_count_90,
    # consec_count_above_95 = consec_count_95
  )
  
  result
}
#Heat
#/*
df_preN1 <- above_thresh(Exp_preN_list1); summary(df_preN1)
df_preN2 <- above_thresh(Exp_preN_list2); summary(df_preN2)
df_preN3 <- above_thresh(Exp_preN_list3); summary(df_preN3)
#*/

df_postN1 <- above_thresh(Exp_postN_list1); summary(df_postN1)
df_postN2 <- above_thresh(Exp_postN_list2); summary(df_postN2)
df_postN3 <- above_thresh(Exp_postN_list3); summary(df_postN3)

#Cold (TDC)
#/*
cold_preN1 <- below_thresh(Exp_preN_list1); summary(cold_preN1)
cold_preN2 <- below_thresh(Exp_preN_list2); summary(cold_preN2)
cold_preN3 <- below_thresh(Exp_preN_list3); summary(cold_preN3)
#*/

cold_postN1 <- below_thresh(Exp_postN_list1); summary(cold_postN1)
cold_postN2 <- below_thresh(Exp_postN_list2); summary(cold_postN2)
cold_postN3 <- below_thresh(Exp_postN_list3); summary(cold_postN3)

#function to calculate mean if garde_alternee
mean_garde_alter <- function(df, name) {
df %>%
  # Create a new "group" column = prefix before underscore (or keep original if no underscore)
  mutate(group = sub("_.*", "", element)) %>%
  group_by(group) %>%
  summarise(
    across(-element, ~ mean(.x, na.rm = TRUE), .names = paste0(name, "_{.col}"))
  ) %>%
    mutate(group=as.numeric(group))
}
#Heat
#/*
df_preN1 <- mean_garde_alter(df_preN1, "preN1")
df_preN2 <- mean_garde_alter(df_preN2, "preN2")
df_preN3 <- mean_garde_alter(df_preN3, "preN3")
#*/

df_postN1 <- mean_garde_alter(df_postN1, "postN1")
df_postN2 <- mean_garde_alter(df_postN2, "postN2")
df_postN3 <- mean_garde_alter(df_postN3, "postN3")

#Cold
#/*
cold_preN1 <- mean_garde_alter(cold_preN1, "preN1")
cold_preN2 <- mean_garde_alter(cold_preN2, "preN2")
cold_preN3 <- mean_garde_alter(cold_preN3, "preN3")
#*/

cold_postN1 <- mean_garde_alter(cold_postN1, "postN1")
cold_postN2 <- mean_garde_alter(cold_postN2, "postN2")
cold_postN3 <- mean_garde_alter(cold_postN3, "postN3")

#/*
dim(df_preN1); dim(df_preN2); dim(df_preN3);  dim(df_postN1); dim(df_postN2); dim(df_postN3); 
dim(cold_preN1); dim(cold_preN2); dim(cold_preN3); dim(cold_postN1); dim(cold_postN2); dim(cold_postN3); 

# save df and cold preN's ----
save(df_preN1, df_preN2, df_preN3,
     cold_preN1, cold_preN2, cold_preN3,
     file=paste0(Exp, "_3monthsPregnancy_ExpMatrices.RData"))
#*/

#length_Exp_preN <- length_Exp_preN[!duplicated(names(length_Exp_preN))]


#POLICY FUNCTION: ONE, TWO, THREE DEGREES EVERYWHERE AND IN EACH PERIOD
#CREATE CF and save ----
#increase of temperature throughout the whole study period, not circumscribed to pregnancy vs. one of the three postnatal periods (unrealistic)
source("Process_data_3monthsPregnancy_BAS_CF.R")
load(paste0(Exp, "_merged_imputed.RData"))
merged_imputed_thisone <- merged_imputed %>% complete(1)
merged_imputed_thisone$id <- merged_imputed$id

Exp_CF1 <- create_db_policy_3monthsPregnancy(id=merged_imputed_thisone$id, Exp, POLICY=1, p90, p95); save(Exp_CF1, file=paste0(Exp,"_Exp_CF1_merged_imputed_3monthsPregnancy.RData"))
Exp_CF2 <- create_db_policy_3monthsPregnancy(id=merged_imputed_thisone$id,Exp, POLICY=2, p90, p95); save(Exp_CF2, file=paste0(Exp,"_Exp_CF2_merged_imputed_3monthsPregnancy.RData"))
Exp_CF3 <- create_db_policy_3monthsPregnancy(id=merged_imputed_thisone$id,Exp, POLICY=3, p90, p95); save(Exp_CF3, file=paste0(Exp,"_Exp_CF3_merged_imputed_3monthsPregnancy.RData"))
