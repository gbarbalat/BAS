create_CoExp <- function(id, CoExp) {
# header ----

load(paste0(CoExp, "_exp.RData"))

#get the CoExp vector for each individual in the final df
all_Exp_preN_list <- get(paste0(CoExp, "_exp_preN_list"))
Exp_preN_list <- unlist(
  lapply(as.character(id), function(i) {
    matches <- grep(paste0("^", i), names(all_Exp_preN_list), value = TRUE)
    all_Exp_preN_list[matches]
  }),
  recursive = FALSE
)

all_Exp_postN_list <- get(paste0(CoExp, "_exp_postN_list"))
Exp_postN_list <- unlist(
  lapply(as.character(id), function(i) {
    matches <- grep(paste0("^", i), names(all_Exp_postN_list), value = TRUE)
    all_Exp_postN_list[matches]
  }),
  recursive = FALSE
)


# preproc step re-CoExp ----


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


## Calculate mean per period ----
calculate_mean <- function(Exp) {
  if (is.null(Exp) || length(Exp) < 10) return(NA_real_)
  mean(Exp %>% unlist, na.rm = TRUE)
}

df_preN <- sapply(Exp_preN_list, calculate_mean); summary(df_preN);
df_postN1 <- sapply(Exp_postN_list1, calculate_mean); summary(df_postN1)
df_postN2 <- sapply(Exp_postN_list2, calculate_mean); summary(df_postN2)
df_postN3 <- sapply(Exp_postN_list3, calculate_mean); summary(df_postN3)


#function to calculate mean if garde_alternee
mean_garde_alter <- function(df, name) {
  
  df <- data.frame(element=names(df), value=df)
    # Create a new "group" column = prefix before underscore (or keep original if no underscore)
  df_summary <- df %>%
    mutate(group = sub("_.*", "", element)) %>%
    group_by(group) %>%
    # summarise(!!paste0(name, "_mean_value") := mean(value, na.rm = TRUE))
    summarise(!!paste0(name, "_", CoExp) := mean(value, na.rm = TRUE))
}

df_preN <- mean_garde_alter(df_preN, "preN") %>% as.data.frame()
df_postN1 <- mean_garde_alter(df_postN1, "postN1") %>% as.data.frame()
df_postN2 <- mean_garde_alter(df_postN2, "postN2") %>% as.data.frame()
df_postN3 <- mean_garde_alter(df_postN3, "postN3") %>% as.data.frame()

all_df <- df_preN %>%
  left_join(df_postN1, by="group") %>%
  left_join(df_postN2, by="group") %>%
  left_join(df_postN3, by="group")
  
return(all_df)

}