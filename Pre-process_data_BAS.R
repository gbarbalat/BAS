# Pre-process aims to create 3 df: Exp (exposure), Out (outcome), Cv (covariates)
# Process aims to merge, filter, select, recode, impute and run perform pre-analysis
# Anal aims to run main, additional and sensitivity analysis

rm(list=ls())
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mice)
library(data.table)
library(stringr)

# header ----
path_to_sas <- "C:/Users/Guillaume/Desktop/PhD_epidemio/Epi/20231130Dem806_585_GB/"
file_ext <- ".sas7bdat"
#Master file
file <- "data_dem806_585_gb"
all_data <- read_sas(paste0(path_to_sas,file,file_ext))


# Exposures ----
exposure_prefix <- c("NO2_preN","NO2_postN", "PM2.5_preN", "PM2.5_postN","PM10_preN", "PM10_postN", 
                     "Tmax_preN", "Tmax_postN","Tmin_preN", "Tmin_postN", "Tmean_preN", "Tmean_postN"
)
# conception = _0 + 365 days
# gestation expressed as WA from M00X_AGEGESTS (AGEGESTS)
# beware conception = 2 weeks of amenorrhea
# previous analyses used a variable called stop_record (e.g. 32 WA)
# this one will use all of the duration of pregnancy and postnatal life till the test is done

#calculate day-wise gestational age
#beware: NA in M00X_AGEGESTJ = 0 day
GESTAGF <- all_data %>%
    #filter(id_Dem806_585_GB %in% Outcome_df$id_Dem806_585_GB) %>%
    select(id_Dem806_585_GB,M00X_AGEGESTS,M00X_AGEGESTJ) %>%
    mutate(M00X_AGEGESTJ=case_when(!is.na(M00X_AGEGESTS) &  is.na(M00X_AGEGESTJ) ~ 
                                     tidyr::replace_na(M00X_AGEGESTJ,0),
                                     .default = as.numeric(M00X_AGEGESTJ)
                                   )
           ) %>%
    mutate(M00X_AGEGEST=M00X_AGEGESTS*7+M00X_AGEGESTJ) 


file_exp1 <- "no2_jour" #from Real et al. coarse resolution
file_exp1 <- "ELFE_BarbalatNO2_NO2_mat_journ_IDGB_20240305.rds" #from Barbalat et al. better resolution
# tmp=readRDS(paste0(path_to_sas,"ELFE_Barbalat_NO2_expo_prenat_IDGB_20240305.rds"))
# tmp=readRDS(paste0(path_to_sas,"ELFE_BarbalatNO2_dates_mat_journ_IDGB_20240305.rds"))
# tmp=readRDS(paste0(path_to_sas,"ELFE_BarbalatNO2_NO2_mat_journ_IDGB_20240305.rds"))

file_exp2 <- "pm10_jour"
file_exp3 <- "pm25_jour"
file_exp4 <- "temp_tmax_jour"
file_exp5 <- "temp_tmin_jour"
file_exp6 <- "temp_tmean_jour"

## NO2 ----
NO2_exp <- readRDS(paste0(path_to_sas,file_exp1))  
#make colnames compatible with other exposure matrices
colnames(NO2_exp)[-1] <- paste0("_", colnames(NO2_exp)[-1])
colnames(NO2_exp)[1] <- "id_Dem806_585_GB"

# garde alternee
idx_mp <- "145373202" # "145373202" "288196197"
NO2_exp %>% filter(id_Dem806_585_GB==paste0(idx_mp,"_p")) %>%
  select(c("_0","_1000","_1500"))
NO2_exp %>% filter(id_Dem806_585_GB==paste0(idx_mp,"_m")) %>%
  select(c("_0","_1000","_1500"))

#add on garde_alter column
p_ <- NO2_exp %>% 
  select(id_Dem806_585_GB) %>%
  filter(grepl("_p$", id_Dem806_585_GB)) %>%
  mutate(garde_alter="_p")

m_ <- NO2_exp %>% 
  select(id_Dem806_585_GB) %>%
  filter(grepl("_m$", id_Dem806_585_GB)) %>%
  mutate(garde_alter="_m")

#check with 295743459 (garde_alter)
NO2_exp <- NO2_exp %>% 
  left_join(p_, by="id_Dem806_585_GB") %>%
  left_join(m_,by="id_Dem806_585_GB") %>%
  mutate(garde_alter=coalesce(garde_alter.x, garde_alter.y)) %>% #garde_alter.x and y never overlap!
  select(-c(garde_alter.x,garde_alter.y)) %>%
  mutate(id_Dem806_585_GB=gsub("_m","",id_Dem806_585_GB)) %>%
  mutate(id_Dem806_585_GB=gsub("_p","",id_Dem806_585_GB)) %>%
  dplyr::mutate(garde_alter = replace_na(garde_alter, "")) %>%
  ####
  #filter(garde_alter!="_p") %>% #take out Dad observations
  ####
  mutate(id_Dem806_585_GB=as.numeric(id_Dem806_585_GB)) %>%
  #join with gestational age dataframe
  left_join(GESTAGF, by = "id_Dem806_585_GB") %>%
  #join with age at test dataframe
  left_join(select(all_data,c(id_Dem806_585_GB,A03F_AGE3A)))
length(unique(NO2_exp$id_Dem806_585_GB)); sum(is.na(NO2_exp$id_Dem806_585_GB)); dim(NO2_exp)
colnames(NO2_exp)[grepl("\\.x$|\\.y$", colnames(NO2_exp))]

#Complete Prenatal exposure 
# Get column names of exposures, excluding ID and non-exposure columns
exposure_cols <- grep("^_\\d+$", colnames(NO2_exp), value = TRUE)

# List containing NO2 exposure for each individual in the prenatal period
NO2_exp_preN_list <- lapply(seq_len(nrow(NO2_exp)), function(i) {
  id_row <- NO2_exp[i, ]
  if (is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  gest_end <- id_row$M00X_AGEGEST - 14 #pregnancy ends date of amenorrhea - 2 weeks
  col_start <- 365#conception starts 365 days post _0
  col_end   <- col_start+gest_end 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(NO2_exp_preN_list) <- paste0(NO2_exp$id_Dem806_585_GB, NO2_exp$garde_alter)

# List containing NO2 exposure for each individual in the postnatal period
NO2_exp_postN_list <- lapply(seq_len(nrow(NO2_exp)), function(i) {

  id_row <- NO2_exp[i, ]
  if (is.na(id_row$A03F_AGE3A) | is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  
  #age at test happened  at least 3 years after birth (365+366+365) and is given in month (assuming that 30 days in a month)
  age_at_test <- (id_row$A03F_AGE3A - 36)*30 + 365+366+365 #age at test in months so transform in days
  col_start <- 365 +  id_row$M00X_AGEGEST - 14 + 1 #postnatal life starts 365 days post _0 (conception) + gestational age - 2 weeks
  col_end   <- col_start+age_at_test 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(NO2_exp_postN_list) <- paste0(NO2_exp$id_Dem806_585_GB, NO2_exp$garde_alter)

#save
save(NO2_exp_preN_list,NO2_exp_postN_list,file="NO2_exp.RData")


## PM10 ----
PM10_exp <- read_sas(paste0(path_to_sas,file_exp2,file_ext)) %>% 
  #filter(id_Dem806_585_GB %in% Outcome_df$id_Dem806_585_GB) %>%
  select(-flag_periode_NA_hors_FRM) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"")) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"NA")) %>%
  mutate_at(vars(all_of(starts_with("_"))), as.numeric)  %>%
  #join with gestational age dataframe
  left_join(GESTAGF, by = "id_Dem806_585_GB") %>%
  #join with age at test dataframe
  left_join(select(all_data,c(id_Dem806_585_GB,A03F_AGE3A)))

idx_mp <- "145373202"  #"145373202" "288196197"
PM10_exp %>% filter(id_Dem806_585_GB==idx_mp) %>%
  select(c("_0","_1000","_1500",garde_alter,id_Dem806_585_GB))
length(unique(PM10_exp$id_Dem806_585_GB)); sum(is.na(PM10_exp$id_Dem806_585_GB)); dim(PM10_exp)
colnames(PM10_exp)[grepl("\\.x$|\\.y$", colnames(PM10_exp))]

#Complete Prenatal exposure 
# Get column names of exposures, excluding ID and non-exposure columns
exposure_cols <- grep("^_\\d+$", colnames(PM10_exp), value = TRUE)

# List containing exposure for each individual in the prenatal period
PM10_exp_preN_list <- lapply(seq_len(nrow(PM10_exp)), function(i) {
  id_row <- PM10_exp[i, ]
  if (is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  gest_end <- id_row$M00X_AGEGEST - 14 #pregnancy ends date of amenorrhea - 2 weeks
  col_start <- 365#conception starts 365 days post _0
  col_end   <- col_start+gest_end 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(PM10_exp_preN_list) <- paste0(PM10_exp$id_Dem806_585_GB, PM10_exp$garde_alter)

# List containing exposure for each individual in the postnatal period
PM10_exp_postN_list <- lapply(seq_len(nrow(PM10_exp)), function(i) {

  id_row <- PM10_exp[i, ]
  if (is.na(id_row$A03F_AGE3A) | is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  
  #age at test happened  at least 3 years after birth (365+366+365) and is given in month (assuming that 30 days in a month)
  age_at_test <- (id_row$A03F_AGE3A - 36)*30 + 365+366+365 #age at test in months so transform in days
  col_start <- 365 +  id_row$M00X_AGEGEST - 14 + 1 #postnatal life starts 365 days post _0 (conception) + gestational age - 2 weeks
  col_end   <- col_start+age_at_test 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(PM10_exp_postN_list) <- paste0(PM10_exp$id_Dem806_585_GB, PM10_exp$garde_alter)

#save
save(PM10_exp_preN_list,PM10_exp_postN_list,file="PM10_exp.RData")


## PM2.5 ----
PM2.5_exp <- read_sas(paste0(path_to_sas,file_exp3,file_ext)) %>%
  #filter(id_Dem806_585_GB %in% Outcome_df$id_Dem806_585_GB) %>%
  select(-flag_periode_NA_hors_FRM) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"")) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"NA")) %>%
  mutate_at(vars(all_of(starts_with("_"))), as.numeric) %>%
  #join with gestational age dataframe
  left_join(GESTAGF, by = "id_Dem806_585_GB") %>%
  #join with age at test dataframe
  left_join(select(all_data,c(id_Dem806_585_GB,A03F_AGE3A)))

colnames(PM2.5_exp)[grepl("\\.x$|\\.y$", colnames(PM2.5_exp))]

#Complete Prenatal exposure 
# Get column names of exposures, excluding ID and non-exposure columns
exposure_cols <- grep("^_\\d+$", colnames(PM2.5_exp), value = TRUE)

# List containing exposure for each individual in the prenatal period
PM2.5_exp_preN_list <- lapply(seq_len(nrow(PM2.5_exp)), function(i) {
  id_row <- PM2.5_exp[i, ]
  if (is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  gest_end <- id_row$M00X_AGEGEST - 14 #pregnancy ends date of amenorrhea - 2 weeks
  col_start <- 365#conception starts 365 days post _0
  col_end   <- col_start+gest_end 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(PM2.5_exp_preN_list) <- paste0(PM2.5_exp$id_Dem806_585_GB, PM2.5_exp$garde_alter)

# List containing exposure for each individual in the postnatal period
PM2.5_exp_postN_list <- lapply(seq_len(nrow(PM2.5_exp)), function(i) {
  
  id_row <- PM2.5_exp[i, ]
  if (is.na(id_row$A03F_AGE3A) | is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  
  #age at test happened  at least 3 years after birth (365+366+365) and is given in month (assuming that 30 days in a month)
  age_at_test <- (id_row$A03F_AGE3A - 36)*30 + 365+366+365 #age at test in months so transform in days
  col_start <- 365 +  id_row$M00X_AGEGEST - 14 + 1 #postnatal life starts 365 days post _0 (conception) + gestational age - 2 weeks
  col_end   <- col_start+age_at_test 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(PM2.5_exp_postN_list) <- paste0(PM2.5_exp$id_Dem806_585_GB, PM2.5_exp$garde_alter)

#save
save(PM2.5_exp_preN_list,PM2.5_exp_postN_list,file="PM2.5_exp.RData")



## Tmax ----
Tmax_exp <- read_sas(paste0(path_to_sas,file_exp4,file_ext)) %>% 
  #filter(id_Dem806_585_GB %in% Outcome_df$id_Dem806_585_GB) %>%
  select(-flag_periode_NA_hors_FRM) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"")) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"NA")) %>%
  mutate_at(vars(all_of(starts_with("_"))), as.numeric) %>%
  #join with gestational age dataframe
  left_join(GESTAGF, by = "id_Dem806_585_GB") %>%
  #join with age at test dataframe
  left_join(select(all_data,c(id_Dem806_585_GB,A03F_AGE3A)))

colnames(Tmax_exp)[grepl("\\.x$|\\.y$", colnames(Tmax_exp))]

#Complete Prenatal exposure 
# Get column names of exposures, excluding ID and non-exposure columns
exposure_cols <- grep("^_\\d+$", colnames(Tmax_exp), value = TRUE)

# List containing exposure for each individual in the prenatal period
Tmax_exp_preN_list <- lapply(seq_len(nrow(Tmax_exp)), function(i) {
  id_row <- Tmax_exp[i, ]
  if (is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  gest_end <- id_row$M00X_AGEGEST - 14 #pregnancy ends date of amenorrhea - 2 weeks
  col_start <- 365#conception starts 365 days post _0
  col_end   <- col_start+gest_end 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(Tmax_exp_preN_list) <- paste0(Tmax_exp$id_Dem806_585_GB, Tmax_exp$garde_alter)

# List containing exposure for each individual in the postnatal period
Tmax_exp_postN_list <- lapply(seq_len(nrow(Tmax_exp)), function(i) {
  
  id_row <- Tmax_exp[i, ]
  if (is.na(id_row$A03F_AGE3A) | is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  
  #age at test happened  at least 3 years after birth (365+366+365) and is given in month (assuming that 30 days in a month)
  age_at_test <- (id_row$A03F_AGE3A - 36)*30 + 365+366+365 #age at test in months so transform in days
  col_start <- 365 +  id_row$M00X_AGEGEST - 14 + 1 #postnatal life starts 365 days post _0 (conception) + gestational age - 2 weeks
  col_end   <- col_start+age_at_test 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(Tmax_exp_postN_list) <- paste0(Tmax_exp$id_Dem806_585_GB, Tmax_exp$garde_alter)

#save
save(Tmax_exp_preN_list,Tmax_exp_postN_list,file="Tmax_exp.RData")



## Tmin ----
Tmin_exp <- read_sas(paste0(path_to_sas,file_exp5,file_ext)) %>% 
  #filter(id_Dem806_585_GB %in% Outcome_df$id_Dem806_585_GB) %>%
  select(-flag_periode_NA_hors_FRM) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"")) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"NA")) %>%
  mutate_at(vars(all_of(starts_with("_"))), as.numeric) %>%
  #join with gestational age dataframe
  left_join(GESTAGF, by = "id_Dem806_585_GB") %>%
  #join with age at test dataframe
  left_join(select(all_data,c(id_Dem806_585_GB,A03F_AGE3A)))

colnames(Tmin_exp)[grepl("\\.x$|\\.y$", colnames(Tmin_exp))]

#Complete Prenatal exposure 
# Get column names of exposures, excluding ID and non-exposure columns
exposure_cols <- grep("^_\\d+$", colnames(Tmin_exp), value = TRUE)

# List containing exposure for each individual in the prenatal period
Tmin_exp_preN_list <- lapply(seq_len(nrow(Tmin_exp)), function(i) {
  id_row <- Tmin_exp[i, ]
  if (is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  gest_end <- id_row$M00X_AGEGEST - 14 #pregnancy ends date of amenorrhea - 2 weeks
  col_start <- 365#conception starts 365 days post _0
  col_end   <- col_start+gest_end 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(Tmin_exp_preN_list) <- paste0(Tmin_exp$id_Dem806_585_GB, Tmin_exp$garde_alter)

# List containing exposure for each individual in the postnatal period
Tmin_exp_postN_list <- lapply(seq_len(nrow(Tmin_exp)), function(i) {
  
  id_row <- Tmin_exp[i, ]
  if (is.na(id_row$A03F_AGE3A) | is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  
  #age at test happened  at least 3 years after birth (365+366+365) and is given in month (assuming that 30 days in a month)
  age_at_test <- (id_row$A03F_AGE3A - 36)*30 + 365+366+365 #age at test in months so transform in days
  col_start <- 365 +  id_row$M00X_AGEGEST - 14 + 1 #postnatal life starts 365 days post _0 (conception) + gestational age - 2 weeks
  col_end   <- col_start+age_at_test 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(Tmin_exp_postN_list) <- paste0(Tmin_exp$id_Dem806_585_GB, Tmin_exp$garde_alter)

#save
save(Tmin_exp_preN_list,Tmin_exp_postN_list,file="Tmin_exp.RData")


## Tmean ----
Tmean_exp <- read_sas(paste0(path_to_sas,file_exp6,file_ext)) %>%
  #filter(id_Dem806_585_GB %in% Outcome_df$id_Dem806_585_GB) %>%
  select(-flag_periode_NA_hors_FRM) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"")) %>%
  mutate_at(vars(all_of(starts_with("_"))), ~ na_if(.,"NA")) %>%
  mutate_at(vars(all_of(starts_with("_"))), as.numeric) %>%
  #join with gestational age dataframe
  left_join(GESTAGF, by = "id_Dem806_585_GB") %>%
  #join with age at test dataframe
  left_join(select(all_data,c(id_Dem806_585_GB,A03F_AGE3A)))

colnames(Tmean_exp)[grepl("\\.x$|\\.y$", colnames(Tmean_exp))]

#Complete Prenatal exposure 
# Get column names of exposures, excluding ID and non-exposure columns
exposure_cols <- grep("^_\\d+$", colnames(Tmean_exp), value = TRUE)

# List containing exposure for each individual in the prenatal period
Tmean_exp_preN_list <- lapply(seq_len(nrow(Tmean_exp)), function(i) {
  id_row <- Tmean_exp[i, ]
  if (is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  gest_end <- id_row$M00X_AGEGEST - 14 #pregnancy ends date of amenorrhea - 2 weeks
  col_start <- 365#conception starts 365 days post _0
  col_end   <- col_start+gest_end 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(Tmean_exp_preN_list) <- paste0(Tmean_exp$id_Dem806_585_GB, Tmean_exp$garde_alter)

# List containing exposure for each individual in the postnatal period
Tmean_exp_postN_list <- lapply(seq_len(nrow(Tmean_exp)), function(i) {
  
  id_row <- Tmean_exp[i, ]
  if (is.na(id_row$A03F_AGE3A) | is.na(id_row$M00X_AGEGEST)) {
    return(NULL)  # or NA, "" etc.
  }
  
  #age at test happened  at least 3 years after birth (365+366+365) and is given in month (assuming that 30 days in a month)
  age_at_test <- (id_row$A03F_AGE3A - 36)*30 + 365+366+365 #age at test in months so transform in days
  col_start <- 365 +  id_row$M00X_AGEGEST - 14 + 1 #postnatal life starts 365 days post _0 (conception) + gestational age - 2 weeks
  col_end   <- col_start+age_at_test 
  # Columns for the individual's gestational window
  exp_cols_this <- paste0("_", col_start:col_end)
  # Keep only columns present in the data
  exp_cols_this <- exp_cols_this[exp_cols_this %in% exposure_cols]
  # Extract values 
  id_row %>%
    select(all_of(exp_cols_this)) 
})
#name the list by IDs
names(Tmean_exp_postN_list) <- paste0(Tmean_exp$id_Dem806_585_GB, Tmean_exp$garde_alter)

#save
save(Tmean_exp_preN_list,Tmean_exp_postN_list,file="Tmean_exp.RData")



# Outcome:  BAS  ----
# A03F_AGE3A
# A03F_SCOREBASB1
# A03F_SCOREBASR0
# A03F_SCOREBASR8
# A03F_SCOREBASR9
# A03F_SCOREBASRAW
# A03F_SCOREBASRAWSUM
# A03F_SCOREBASABIL
# A03F_SCOREBASPERC
# A03F_APREBAS

Out <- all_data %>%
  select(id_Dem806_585_GB
         , A03F_SCOREBASB1#reponses correctes
         , A03F_SCOREBASR0#rep incorrectes
         , A03F_SCOREBASR8#refus
         , A03F_SCOREBASR9#nsp
         , A03F_SCOREBASRAW
         , A03F_SCOREBASRAWSUM#maximum theorique de planches
         , A03F_SCOREBASABIL
         , A03F_SCOREBASPERC#pct du score d'aptitude
         , A03F_APREBAS
  ) %>%
  mutate(ALL=A03F_SCOREBASRAW)
length(unique(Out$id_Dem806_585_GB)); sum(is.na(Out$id_Dem806_585_GB)); dim(Out)
save(Out, file = "Out.RData")


# CoVar from add. files ----

## Moves ----
file_Moves <- "DEMENAG"
Moves <- read.csv(paste0(path_to_sas,file_Moves,".csv"))
str(Moves); 
table(Moves$A01M_DEMENAG, useNA = "always"); table(Moves$A01P_DEMENAG, useNA = "always");table(Moves$A01R_DEMENAG, useNA = "always");
table(Moves$A02M_DEMENAG, useNA = "always"); table(Moves$A02P_DEMENAG, useNA = "always");table(Moves$A02R_DEMENAG, useNA = "always");
table(Moves$A03R_DEMENAG, useNA = "always");
Moves <- Moves %>%
  mutate(
    A01R_DEMENAG = if_else(is.na(A01M_DEMENAG), A01P_DEMENAG, A01M_DEMENAG),
    A02R_DEMENAG = if_else(is.na(A02M_DEMENAG), A02P_DEMENAG, A02M_DEMENAG)
  )
Moves <- Moves %>%
  select(id_Dem806_585_GB, A01R_DEMENAG, A02R_DEMENAG, A03R_DEMENAG)

## EDI ----
file_EDI2007 <- "edi_score_jour_2007"
file_EDI2011 <- "edi_score_jour_2011"
file_EDI2015 <- "edi_score_jour_2015"

EDI_2007 <- read_sas(paste0(path_to_sas,file_EDI2007,file_ext)) %>% #select(-garde_alter) %>%
  filter(garde_alter!="_p") %>%
  select(c(id_Dem806_585_GB,"_0")) %>% rename(EDI="_0") %>% 
  mutate(EDI=na_if(EDI,"")) %>%
  mutate(EDI=na_if(EDI,"NA")) %>%
  mutate_all(as.numeric) %>%
  rename(EDI_2007=EDI)

EDI_2011 <- read_sas(paste0(path_to_sas,file_EDI2011,file_ext)) %>% #select(-garde_alter) %>%
  filter(garde_alter!="_p") %>%
  select(c(id_Dem806_585_GB,"_0")) %>% rename(EDI="_0") %>% 
  mutate(EDI=na_if(EDI,"")) %>%
  mutate(EDI=na_if(EDI,"NA")) %>%
  mutate_all(as.numeric)  %>%
  rename(EDI_2011=EDI)

EDI_2015 <- read_sas(paste0(path_to_sas,file_EDI2015,file_ext)) %>% #select(-garde_alter) %>%
  filter(garde_alter!="_p") %>%
  select(c(id_Dem806_585_GB,"_0")) %>% rename(EDI="_0") %>% 
  mutate(EDI=na_if(EDI,"")) %>%
  mutate(EDI=na_if(EDI,"NA")) %>%
  mutate_all(as.numeric)  %>%
  rename(EDI_2015=EDI)
length(unique(EDI_2015$id_Dem806_585_GB)); sum(is.na(EDI_2015$id_Dem806_585_GB)); dim(EDI_2015)


## NDVI ----
file_NDVI <- "ndvi"
NDVI_var <- read_sas(paste0(path_to_sas,file_NDVI,file_ext)) %>% #select(-garde_alter) %>%
  filter(year %in% c(2010,2011,2012,2013,2014),
         garde_alter!="_p") %>%
  mutate(buffer100m_pctdata_ete=na_if(buffer100m_pctdata_ete,"")) %>%
  mutate(buffer100m_pctdata_ete=na_if(buffer100m_pctdata_ete,"NA")) %>%
  mutate(buffer100m_pctdata_ete=as.numeric(buffer100m_pctdata_ete)) %>%
  mutate(buffer100m_ndvi_ete=case_when(buffer100m_pctdata_ete < 0.75 ~ NA_character_, 
                                       TRUE ~ buffer100m_ndvi_ete)) %>%
  
  select(c(id_Dem806_585_GB, year, buffer100m_ndvi_ete)) %>%
  mutate(buffer100m_ndvi_ete=na_if(buffer100m_ndvi_ete,"")) %>%
  mutate(buffer100m_ndvi_ete=na_if(buffer100m_ndvi_ete,"NA")) %>%
  
  mutate(buffer100m_ndvi_ete=as.numeric(buffer100m_ndvi_ete)) %>%

  #summarise when for a single subject indicates same year multiple times
  group_by(id_Dem806_585_GB, year) %>%
  summarise(buffer100m_ndvi_ete = mean(buffer100m_ndvi_ete, na.rm = TRUE)) %>%
  ungroup() %>%
  
  #make it wide
  pivot_wider(
    id_cols = id_Dem806_585_GB,
    names_from = year,
    values_from = buffer100m_ndvi_ete,
    names_prefix = "buffer100m_ndvi_ete_"
  )
length(unique(NDVI_var$id_Dem806_585_GB)); sum(is.na(NDVI_var$id_Dem806_585_GB)); dim(NDVI_var)


## EQR12 ----
#eqr1 #Alimentation de la mère pendant les 3 derniers mois de la grossesse (EQR1)
#eqr3 #Alimentation lactée du nourrisson (EQR3)
#eqr12 #Socio-demo
file_eqr12 <- "eqr12"
EQR12 <- read_sas(paste0(path_to_sas,file_eqr12,file_ext)) %>%
  select(c(id_Dem806_585_GB,
           meduc_2m,feduc_2m,
           revenu_part_qui_2m,revenu_part_qui_1y,revenu_part_qui_2y,revenu_part_qui_3y,
           
           "c_emp_2m","c_emp_1y","c_emp_2y","c_emp_3y", #ceux qui cohabitent ont une activite pro #"M00M2_CSP1M", "M00M2_CSP1P",
           "sib_2m","sib_1y","sib_2y","sib_3y", #siblings
           "Child_hhld_2m","child_hhld_1y","child_hhld_2y","child_hhld_3y",#Ou vit l'enfant #"A02M_EFVIT",
           "house_ownership_2m","house_ownership_1y","house_ownership_2y","house_ownership_3y", #
           
           #only one measure
           "mimm", "fimm"# migration status # "M00M2_LIEUNAISM",# born in France or another country
                      )
         )
length(unique(EQR12$id_Dem806_585_GB)); sum(is.na(EQR12$id_Dem806_585_GB)); dim(EQR12)


##eqr29 ----
#eqr29, #Langues parlees ---
# f) La variables NBLANGMEN donne le nombre de langues différentes parlées
# dans le ménage (synthèse des langues parlées à l’enfant par chacun des parents
# et des autres langues parlées dans le ménage déclarées par chaque parent)
file_eqr29 <- "eqr29"
NBANGMEN <- read_sas(paste0(path_to_sas,file_eqr29,file_ext)) %>% 
  select(c(id_Dem806_585_GB,starts_with("NBLANGMEN_"))) %>%
  select(c(id_Dem806_585_GB,NBLANGMEN_2m)) 
length(unique(NBANGMEN$id_Dem806_585_GB)); sum(is.na(NBANGMEN$id_Dem806_585_GB)); dim(NBANGMEN)


## Urban ----
#bc4_mater_mere
#bc4_2mois_mere
#bc4_1an_mere
file_bc4 <- "bc4_1an_mere" #bc4_1an_mere
bc4_mater <- read_sas(paste0(path_to_sas,"bc4_1an_mere",file_ext)) %>% #%>% select(-garde_alter) %>%
  select(-M00X4_CATAEUR2010) %>%
  mutate(M00X4_TAU2010=as.character(M00X4_TAU2010)) 
length(unique(bc4_mater$id_Dem806_585_GB)); sum(is.na(bc4_mater$id_Dem806_585_GB)); dim(bc4_mater)
#TAU2010  taille de l’aire urbaine
#voir doc BdDC_VariablesUrbanisation_2011


#Merge all covar
Cv <- all_data %>%
  select(-c(A03F_SCOREBASB1#reponses correctes
              , A03F_SCOREBASR0#rep incorrectes
              , A03F_SCOREBASR8#refus
              , A03F_SCOREBASR9#nsp
              , A03F_SCOREBASRAW
              , A03F_SCOREBASRAWSUM#maximum theorique de planches
              , A03F_SCOREBASABIL
              , A03F_SCOREBASPERC#pct du score d'aptitude
              , A03F_APREBAS)
         ) %>% 
  left_join(Moves,by="id_Dem806_585_GB") %>%
  left_join(EQR12,by="id_Dem806_585_GB") %>%
  left_join(EDI_2007,by="id_Dem806_585_GB") %>%
  left_join(EDI_2011,by="id_Dem806_585_GB") %>%
  left_join(EDI_2015,by="id_Dem806_585_GB") %>%
  left_join(NDVI_var,by="id_Dem806_585_GB") %>%
  left_join(NBANGMEN,by="id_Dem806_585_GB") %>%
  left_join(bc4_mater,by="id_Dem806_585_GB")
length(unique(Cv$id_Dem806_585_GB)); sum(is.na(Cv$id_Dem806_585_GB)); dim(Cv)
save(Cv, file = "Cv.RData")
