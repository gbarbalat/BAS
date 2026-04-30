#take merged_imputed and simply change Heat and Cold variables
#Or re-do everything (inc. imputation step)
#in a +2degree world , all other things being equal is not necessarily true .... perhaps there will be some impact on illnesses, birth and pregnancy health outcomes


create_db_policy <- function(id, Exp, POLICY, p90, p95) {
# header ----
#POLICY <- 1 #1, 2 or 3 degrees

path_to_sas <- "C:/Users/Guillaume/Desktop/PhD_epidemio/Epi/20231130Dem806_585_GB/"
file_ext <- ".sas7bdat"

load(paste0(Exp, "_exp.RData"))

#get the exposure vector for each individual in the final df
all_Exp_preN_list <- get(paste0(Exp, "_exp_preN_list"))
Exp_preN_list <- unlist(
  lapply(as.character(id), function(i) {
    matches <- grep(paste0("^", i), names(all_Exp_preN_list), value = TRUE)
    all_Exp_preN_list[matches]
  }),
  recursive = FALSE
)

all_Exp_postN_list <- get(paste0(Exp, "_exp_postN_list"))
Exp_postN_list <- unlist(
  lapply(as.character(id), function(i) {
    matches <- grep(paste0("^", i), names(all_Exp_postN_list), value = TRUE)
    all_Exp_postN_list[matches]
  }),
  recursive = FALSE
)

# apply policy ----
Exp_preN_list <- lapply(Exp_preN_list, function(x) x + POLICY)
Exp_postN_list <- lapply(Exp_postN_list, function(x) x + POLICY)

# preproc step re- exposure ----


## split postN into 3 periods ----
# Function to extract a period
extract_period <- function(x, start, end) {
  n <- length(x)
  if (is.null(x) | start>n) return(NULL)  # no data in this range
  x[start:min(end, n)]
}

# Create new lists
Exp_postN_list1 <- lapply(Exp_postN_list, extract_period, start = 1, end = 365)
Exp_postN_list2 <- lapply(Exp_postN_list, extract_period, start = 366, end = 365+366)
Exp_postN_list3 <- lapply(Exp_postN_list, extract_period, start = 365+366+1, end = Inf)

# Preserve names
names(Exp_postN_list1) <- names(Exp_postN_list)
names(Exp_postN_list2) <- names(Exp_postN_list)
names(Exp_postN_list3) <- names(Exp_postN_list)


## Ndays above thresh ----
above_thresh <- function(Exp) {
  
  # Function to count temperature above threshold
  count_above <- function(x, threshold) {
    if (is.null(x) | length(x) < 10) return(NA)  # no data in this range
    
    sum(x > threshold, na.rm = TRUE)
  }
  
  # Count consecutive days above threshold
  count_above_consecutive <- function(x, threshold) {
    if (is.null(x) |  length(x) < 10) return(NA)  # no data in this range
    
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

#Heat
df_preN <- above_thresh(Exp_preN_list); summary(df_preN); 
df_postN1 <- above_thresh(Exp_postN_list1); summary(df_postN1)
df_postN2 <- above_thresh(Exp_postN_list2); summary(df_postN2)
df_postN3 <- above_thresh(Exp_postN_list3); summary(df_postN3)


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
df_preN <- mean_garde_alter(df_preN, "preN")
df_postN1 <- mean_garde_alter(df_postN1, "postN1")
df_postN2 <- mean_garde_alter(df_postN2, "postN2")
df_postN3 <- mean_garde_alter(df_postN3, "postN3")

return(cbind(df_preN,df_postN1,df_postN2,df_postN3))

}