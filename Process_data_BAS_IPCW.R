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
args<- 1

#mice parameters
m <- 10
maxit <- 20
m <- 2
maxit <- 1

path_to_sas <- "/bettik/barbalag/BAS/"
path_to_sas <- "C:/Users/Guillaume/Desktop/PhD_epidemio/Epi/20231130Dem806_585_GB/"

file_ext <- ".sas7bdat"
all_Exp <- c("Tmin", "Tmax", "Tmean"); Exp <- all_Exp[as.numeric(args[1])]

load("Out.RData")
load("Cv.RData")
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
Exp_postN_list1 <- lapply(Exp_postN_list, extract_period, start = 1, end = 365)
Exp_postN_list2 <- lapply(Exp_postN_list, extract_period, start = 366, end = 365+366)
Exp_postN_list3 <- lapply(Exp_postN_list, extract_period, start = 365+366+1, end = Inf)

# Preserve names
names(Exp_postN_list1) <- names(Exp_postN_list)
names(Exp_postN_list2) <- names(Exp_postN_list)
names(Exp_postN_list3) <- names(Exp_postN_list)

## Length of all elements in Exp ----
length_Exp_preN <- sapply(Exp_preN_list, length)
length_Exp_postN3 <- sapply(Exp_postN_list, length) - (365+366)
#save(length_Exp_preN, length_Exp_postN3, file=paste0(Exp, "_length_Exp.RData"))


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
df_preN <- above_thresh(Exp_preN_list); summary(df_preN); 
df_postN1 <- above_thresh(Exp_postN_list1); summary(df_postN1)
df_postN2 <- above_thresh(Exp_postN_list2); summary(df_postN2)
df_postN3 <- above_thresh(Exp_postN_list3); summary(df_postN3)

#Cold (TDC)
cold_preN <- below_thresh(Exp_preN_list); summary(cold_preN); 
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
df_preN <- mean_garde_alter(df_preN, "preN")
df_postN1 <- mean_garde_alter(df_postN1, "postN1")
df_postN2 <- mean_garde_alter(df_postN2, "postN2")
df_postN3 <- mean_garde_alter(df_postN3, "postN3")

#Cold
cold_preN <- mean_garde_alter(cold_preN, "preN")
cold_postN1 <- mean_garde_alter(cold_postN1, "postN1")
cold_postN2 <- mean_garde_alter(cold_postN2, "postN2")
cold_postN3 <- mean_garde_alter(cold_postN3, "postN3")

# merged_ ----
# merge Exp, Out and Cv.  
# explore #, uniqueness and NA in ID; colnames with .x or .y
dim(Out); dim(Cv); 
dim(df_preN); dim(df_postN1); dim(df_postN2); dim(df_postN3); 
dim(cold_preN); dim(cold_postN1); dim(cold_postN2); dim(cold_postN3); 

merged_ <- Out %>%
  
  left_join(df_preN, by = c("id_Dem806_585_GB" = "group")) %>%
  left_join(df_postN1, by = c("id_Dem806_585_GB" = "group")) %>%
  left_join(df_postN2, by = c("id_Dem806_585_GB" = "group")) %>%
  left_join(df_postN3, by = c("id_Dem806_585_GB" = "group")) %>%
  
  left_join(cold_preN, by = c("id_Dem806_585_GB" = "group")) %>%
  left_join(cold_postN1, by = c("id_Dem806_585_GB" = "group")) %>%
  left_join(cold_postN2, by = c("id_Dem806_585_GB" = "group")) %>%
  left_join(cold_postN3, by = c("id_Dem806_585_GB" = "group")) %>%
  
  left_join(Cv, by="id_Dem806_585_GB")
length(unique(merged_$id_Dem806_585_GB)); sum(is.na(merged_$id_Dem806_585_GB)); dim(merged_)
colnames(merged_)[grepl("\\.x$|\\.y$", colnames(merged_))]


# merged_col_obs ----
merged_col_obs <- merged_

# obvious select, 
#M00X_HTA No, chronique, previous Pregnancy (PP), recoupe #M00X_HTAG No, with or w/o proteinurie
#M00X_DIABETE 4 categories no, type I, type II, PP, overlaps perfectly with #M00X_DIABGEST (No/Yes)
# risk is high that diag gest / HTA grossesse only are recorded (disorder of pregnancy -- mediator)
variable_set_1 <- c(
  
  #Outcome
  "ALL",
  
  #Exposure
  "preN_count_above_90", "preN_consec_count_above_90", "preN_count_above_95", "preN_consec_count_above_95",
  "postN1_count_above_90", "postN1_consec_count_above_90", "postN1_count_above_95", "postN1_consec_count_above_95",
  "postN2_count_above_90", "postN2_consec_count_above_90", "postN2_count_above_95", "postN2_consec_count_above_95",
  "postN3_count_above_90", "postN3_consec_count_above_90", "postN3_count_above_95", "postN3_consec_count_above_95",
  
  #Cold TDC
  "preN_count_below_10", "preN_consec_count_below_10", "preN_count_below_05", "preN_consec_count_below_05",
  "postN1_count_below_10", "postN1_consec_count_below_10", "postN1_count_below_05", "postN1_consec_count_below_05",
  "postN2_count_below_10", "postN2_consec_count_below_10", "postN2_count_below_05", "postN2_consec_count_below_05",
  "postN3_count_below_10", "postN3_consec_count_below_10", "postN3_count_below_05", "postN3_consec_count_below_05",
  
  # EDI, NDVI, Urban
  "EDI_2007","EDI_2011","EDI_2015",
  "buffer100m_ndvi_ete_2010","buffer100m_ndvi_ete_2011","buffer100m_ndvi_ete_2012","buffer100m_ndvi_ete_2013","buffer100m_ndvi_ete_2014", # 
  "M00X4_TAU2010",
  
  # CSP mum and dad
  "M00M2_CSP1M", 
  "M00M2_CSP1P",
  
  # From EQR12 (Socio-Demo)
  "c_emp_2m","c_emp_1y","c_emp_2y","c_emp_3y", #ceux qui cohabitent ont une activite pro 
  "sib_2m","sib_1y","sib_2y","sib_3y", #siblings
  "Child_hhld_2m","child_hhld_1y","child_hhld_2y","child_hhld_3y",#Ou vit l'enfant #"A02M_EFVIT",
  "revenu_part_qui_2m","revenu_part_qui_1y","revenu_part_qui_2y","revenu_part_qui_3y", # 
  "house_ownership_2m","house_ownership_1y","house_ownership_2y","house_ownership_3y", #
  
  # From EQR49 (Languages)
  "NBLANGMEN_2m", # "NBLANGMEN_1y",  "NBLANGMEN_2y",  "NBLANGMEN_3y", #nb of languages spoken at home
  
  # Moves
  "A01R_DEMENAG","A02R_DEMENAG","A03R_DEMENAG",
  
  #only one measure
  #"educ_2m",  ### #higher of meduc & feduc,
  "feduc_2m","meduc_2m",
  #"imm", # migration status # 
  "mimm","fimm",
  
  "M00M1_VAGUE", #season of birth
  #region
  
  "SEXE_ENF",
  
  "M00M2_AGEM", #mother's age at birth
  "M00M2_AGEP", #father's age at birth
  
  "M00M2_ENFGANT", 
  "M00M2_FAF",

  #"M00M2_BMIMAVTG",
  "M00M2_POIMAVTG",
  "M00M2_TAIM",
  
  "M00M2_FQALCOOL", #freq etOH (ordinal)
  
  "M00M3_FQCAFE", #6 categ
  
  "M00M3_POISGEN",# Never, <1/31, 1-3/31, 1/7, 2-5/7, 7/7, xxx7/7, always
  
  # "M00M3_VITB9", # No Yes DNK
  'M00M3_VITAB9',
  'M00M2_ACIDEFOL',
  
  "M00M3_MGOMEGA3",# Never, <1/7, x/7, ~7/7, Allways
  
  #### 2M
  "M02M_TYPALI", # = breast, breast + bottle, bottle
  
  #tobacco
  #"TOBACCO",
  #"M02_EXPTAB", #   
  "M02M_EXPTAB", # exposure to tobacco 5 categ
  "M02P_EXPTAB", #
  "M00M2_TABAG",#  
  
  # "M00M2_EXPTAB", ##passive tobacco (home) 
  "M00M2_EXPTABD", 
  "M00M2_EXPTABLF", 
  
  #### 1Y
  #"A01M_HxNDD",
  #"A01P_HxNDD",
    'A01M_DIFMATH',
    'A01M_DIFLIR',
    'A01M_DIFORTH',
    'A01M_RLGG',
    'A01M_DIFORA',
    'A01M_PBCOM',
    'A01P_DIFMATH',
    'A01P_DIFLIR',
    'A01P_DIFORTH',
    'A01P_RLGG',
    'A01P_DIFORA',
    'A01P_PBCOM',
  
  #### 3Y
  "A03F_AGE3A"# age in months
  
)
merged_col_obs <- merged_col_obs %>%
  select(all_of(variable_set_1), id_Dem806_585_GB)


# filter (excl criteria) 
# exclude multiple pregnancies
# exclude no outcome
# exclude no exposure
#this is the multiple pregnancy file
file_rangalea <- 'rangalea'
RANGALEA <- read_sas(paste0(path_to_sas,file_rangalea,file_ext))

#id multiple pregn
id_MB <- RANGALEA %>%
  filter(!is.na(M00M1_RANGALEA)) %>%
  select(id_Dem806_585_GB) 
#remain in main df
# merged_col_obs <- merged_col_obs %>%
#   filter(!id_Dem806_585_GB %in% id_MB$id_Dem806_585_GB)

#id_3y_NoQS from remaining participants
id_3y_NoQS <- merged_col_obs %>%
  filter(is.na(ALL)) %>%
  select(id_Dem806_585_GB) 
#remain in main df
# merged_col_obs <- merged_col_obs %>%
#   filter(!id_Dem806_585_GB %in% id_3y_NoQS$id_Dem806_585_GB)

#id NA in exposure from remaining participants

#NULL in Exp will be NA 
null_preN  <- names(Exp_preN_list)[sapply(Exp_preN_list, is.null)]
null_postN  <- names(Exp_postN_list)[sapply(Exp_postN_list, is.null)]

# Identify other problematic elements
has_na <- sapply(Exp_preN_list, function(x) length(x) > 0 && any(is.na(x)))
all_na <- sapply(Exp_preN_list, function(x) length(x) > 0 && all(is.na(x)))

# Count and id
length_at_least_one_NA <- sum(has_na); names_with_na <- names(Exp_preN_list)[has_na]
length_all_NA <- sum(all_na); names_all_na <- names(Exp_preN_list)[all_na]
all_names_NA_preN <- unique(c(names_all_na, names_with_na))

# same with postnatal exposure
# Identify other problematic elements
has_na <- sapply(Exp_postN_list, function(x) length(x) > 0 && any(is.na(x)))
all_na <- sapply(Exp_postN_list, function(x) length(x) > 0 && all(is.na(x)))

# Count and id
length_at_least_one_NA <- sum(has_na); names_with_na <- names(Exp_postN_list)[has_na]
length_all_NA <- sum(all_na); names_all_na <- names(Exp_postN_list)[all_na]
all_names_NA_postN <- unique(c(names_all_na, names_with_na))

# id_NAExp from remaining participants
# concatenate preN and postN
id_NAExp <- sub("_.*", "",c(null_preN, null_postN, all_names_NA_preN, all_names_NA_postN)) %>% unique %>% as.numeric

#remain in main df
# merged_col_obs <- merged_col_obs %>%
#   filter(!id_Dem806_585_GB %in% id_NAExp)

#/*
id_MB;id_3y_NoQS;id_NAExp
#*/

# 1st step of recoding - inc. na_if, make categ, ordinal and num
merged_col_obs <- merged_col_obs %>% 
    
    #Using rowwise() when it's not needed will just slow down your code. It does not affect correctness with vectorized functions like case_when() or 
    #recode(), but it’s unnecessary and inefficient.
    rowwise() %>%
    # education: define a new variable, te max from father and mother - rowwise necessary
    mutate(educ_2m=case_when(is.na(meduc_2m) & is.na(feduc_2m) ~ NA_integer_,
                             is.na(meduc_2m) & !is.na(feduc_2m) ~ feduc_2m,
                             !is.na(meduc_2m) & is.na(feduc_2m) ~ meduc_2m,
                             !is.na(meduc_2m) & !is.na(feduc_2m) ~ max(c(feduc_2m,meduc_2m))
                             ),
           
    #immigration status
    imm=case_when(mimm==1 | fimm==1 ~ "atleastone_imm",
                  mimm==4 & fimm==4 ~ "both_fra_fra",   
                  mimm==2 | mimm==3 | fimm==2 | fimm==3 ~ "both_fra_atleastoneparentimm",
                  TRUE ~ NA_character_
                  ),
    
    #Move
    A01R_DEMENAG=recode(A01R_DEMENAG,
                        "1"="Yes",
                        "2"="No"),
    A02R_DEMENAG=recode(A02R_DEMENAG,
                        "1"="Yes",
                        "2"="No"),
    A03R_DEMENAG=recode(A03R_DEMENAG,
                        "1"="Yes",
                        "2"="No"),
    
    #"c_emp_2m","c_emp_1y","c_emp_2y","c_emp_3y", #ceux qui cohabitent ont une activite pro 
    c_emp_2m=recode(c_emp_2m,
                    "1"="1",
                    "2"="None",
                    "3"="both"),
    c_emp_1y=recode(c_emp_1y,
                    "1"="1",
                    "2"="None",
                    "3"="both"),
    c_emp_2y=recode(c_emp_2y,
                    "1"="1",
                    "2"="None",
                    "3"="both"),
    c_emp_3y=recode(c_emp_3y,
                    "1"="1",
                    "2"="None",
                    "3"="both"),
    
    #"sib_2m","sib_1y","sib_2y","sib_3y", #siblings
    sib_2m=recode(sib_2m,
                  "0"="0",
                  "1"="1",
                  "2"="2",
                  .default = "3+"),
    sib_1y=recode(sib_1y,
                  "0"="0",
                  "1"="1",
                  "2"="2",
                  .default = "3+"),
    sib_2y=recode(sib_2y,
                  "0"="0",
                  "1"="1",
                  "2"="2",
                  .default = "3+"),
    sib_3y=recode(sib_3y,
                  "0"="0",
                  "1"="1",
                  "2"="2",
                  .default = "3+"),
    
    #"Child_hhld_2m","child_hhld_1y","child_hhld_2y","child_hhld_3y",#Ou vit l'enfant #"A02M_EFVIT",
    Child_hhld_2m=recode(Child_hhld_2m,
                  "1"="2parents",
                  .default = "not_2parents"),
    child_hhld_1y=recode(child_hhld_1y,
                         "1"="2parents",
                         .default = "not_2parents"),
    child_hhld_2y=recode(child_hhld_2y,
                         "1"="2parents",
                         .default = "not_2parents"),
    child_hhld_3y=recode(child_hhld_3y,
                         "1"="2parents",
                         .default = "not_2parents"),
    
    #"house_ownership_2m","house_ownership_1y","house_ownership_2y","house_ownership_3y", #
    house_ownership_2m=recode(house_ownership_2m,
                         "1"="Landlord",
                         "2"="Rent",
                         "3"="Free_Rent_HLM",
                         "4"="Free_Rent_HLM",
                         "5"="Rent",
                         .default = NA_character_),
    house_ownership_1y=recode(house_ownership_1y,
                              "1"="Landlord",
                              "2"="Rent",
                              "3"="Free_Rent_HLM",
                              "4"="Free_Rent_HLM",
                              "5"="Rent",
                              .default = NA_character_),
    house_ownership_2y=recode(house_ownership_2y,
                              "1"="Landlord",
                              "2"="Rent",
                              "3"="Free_Rent_HLM",
                              "4"="Free_Rent_HLM",
                              "5"="Rent",
                              .default = NA_character_),
    house_ownership_3y=recode(house_ownership_3y,
                              "1"="Landlord",
                              "2"="Rent",
                              "3"="Free_Rent_HLM",
                              "4"="Free_Rent_HLM",
                              "5"="Rent",
                              .default = NA_character_),
    
    #CSP
    M00M2_CSP1P=recode(M00M2_CSP1P,
                       "1"="Agr_Stud_None",
                       "7"="Agr_Stud_None",
                       "8"="Agr_Stud_None",
                       "9"=NA_character_,
                       .default = as.character(M00M2_CSP1P)
    ),
    M00M2_CSP1M=recode(M00M2_CSP1M,
                       "1"="Agr_Stud_None",
                       "7"="Agr_Stud_None",
                       "8"="Agr_Stud_None",
                       "9"=NA_character_,
                       .default = as.character(M00M2_CSP1P)
                       ),
    
    #exposure to tobacco at 2 months - rowwise necessary 
    M02_EXPTAB=case_when(is.na(M02M_EXPTAB) & is.na(M02P_EXPTAB) ~ NA_integer_,
                             is.na(M02M_EXPTAB) & !is.na(M02P_EXPTAB) ~ M02P_EXPTAB,
                             !is.na(M02M_EXPTAB) & is.na(M02P_EXPTAB) ~ M02M_EXPTAB,
                             !is.na(M02M_EXPTAB) & !is.na(M02P_EXPTAB) ~ max(c(M02M_EXPTAB,M02P_EXPTAB))
    ),
    #Never; <1h/day; 1-2 hrs; 2-5 hrs; >5hrs
    
    #exposure to tobacco during pregnancy
    M00M2_EXPTAB=case_when(M00M2_EXPTABD %in% c(1,2,3,4) | 
                             M00M2_EXPTABLF %in% c(1,2,3,4) ~ "Yes",#passive home or  other than home
                           M00M2_EXPTABD==0 & M00M2_EXPTABLF==0 ~ "No",
                           .default=NA_character_),
    #Never; <1h/day; 1-2 hrs; 2-5 hrs; >5hrs
    
    SEXE_ENF=recode(SEXE_ENF,"1"="Male", "2"= "Female"),
    
    #group if 3 or more languages
    NBLANGMEN_2m=recode(NBLANGMEN_2m,
                        "3"="3_",
                        "4"="3_",
                        "5"="3_",
                        "6"="3_",
                        .default = as.character(NBLANGMEN_2m)
    ),
    
    #if you wanted to categorize mum and dad's age
    # M00M2_AGEM=cut(M00M2_AGEM,breaks=c(0,25,30,35,40,80)),
    # M00M2_AGEP=cut(M00M2_AGEP,breaks=c(0,25,30,35,40,45,80)),
    
    #previous pregnancies: define a new variable from 2 other variables
    M00M2_ENFGANT=case_when(is.na(M00M2_FAF) ~ NA_character_,
                            is.na(M00M2_ENFGANT) ~ "No",
                            M00M2_ENFGANT==0 ~"No",
                            M00M2_ENFGANT==1 ~ "Yes",
                            M00M2_ENFGANT==9 ~ NA_character_),
    
    #BMI: define a new variable 
    M00M2_BMIMAVTG=M00M2_POIMAVTG/(M00M2_TAIM/100)^2,

    #etOH during pregnancy
    M00M2_FQALCOOL=recode(M00M2_FQALCOOL, 
                          "2"="More_1m",
                          "3"="More_1m",
                          "4"="More_1m",
                          "5"="More_1m",
                          "7"=NA_character_,
                          .default=as.character(M00M2_FQALCOOL)),
    #M00M2_FQALCOOL: Never, <1/31, 2-4/31, {2-3/7,>4/7,7/7}, B/f pregn, DNK
    
    # M00M3_FQCAFE=recode(M00M3_FQCAFE,
    #                     "0"="Less_1_D","1"="Less_1_D","2"="Less_1_D","3"="Less_1_D",
    #                     "4"="More_1_D", "5"="More_1_D","6"="More_1_D",
    #                     .default=NA_character_),
    # #0, <1/M, 1-3/M, 1/7, 2-5/7, 7/7, xxx/7
    #"M00M3_CAFE", #Never, More, Same, Less
    #"M00M3_FQCAFE", #6 categ
    
    M00M3_POISGEN=recode(M00M3_POISGEN,
                         "4"="More_2_W", "5"="More_2_W","6"="More_2_W",
                         "7"=NA_character_,
                         .default = as.character(M00M3_POISGEN)),
    #"M00M3_POISGEN",# Never, <1/31, 1-3/31, 1/7, 2-5/7, 7/7, xxx7/7, DNK
    # M00M3_CONSOPOIS More than once a month in the third trimester No, Yes
    
    #VITB9 define a new var from VITB9 and FOLIC ACID variables
    M00M3_VITB9=case_when(M00M3_VITAB9==1 | M00M2_ACIDEFOL==1 ~ "Yes",#M00M3_VITAG==1 | M00M3_VITAMEL==1 | 
                          M00M2_ACIDEFOL==2 ~ "No",#M00M3_VITAG==0 & M00M3_VITAMEL==0 &
                          M00M2_ACIDEFOL==9 ~ NA_character_,
                          .default=NA_character_),
    
    
    M00M3_MGOMEGA3=recode(M00M3_MGOMEGA3, #passive other than home
                          "3"="7_7","4"="7_7",
                          "7"=NA_character_,
                          .default=as.character(M00M3_MGOMEGA3)),
    # "M00M3_ACIDEG", # Yes
    # "M00M3_MGOMEGA3",# Never, <1/7, x/7, ~7/7, Always
    
    # hx of neurodev disorder in mum or dad
    A01M_HxNDD=case_when(A01M_DIFMATH==1 | A01M_DIFLIR==1 | A01M_DIFORTH==1 | A01M_RLGG==1 | A01M_DIFORA==1 | A01M_PBCOM==1 ~ "Yes",
                         A01M_DIFMATH==2 & A01M_DIFLIR==2 & A01M_DIFORTH==2 & A01M_RLGG==2 & A01M_DIFORA==2 & A01M_PBCOM==2 ~ "No",
                         .default = NA_character_),
    
    A01P_HxNDD=case_when(A01P_DIFMATH==1 | A01P_DIFLIR==1 | A01P_DIFORTH==1 | A01P_RLGG==1 | A01P_DIFORA==1 | A01P_PBCOM==1 ~ "Yes",
                         A01P_DIFMATH==2 & A01P_DIFLIR==2 & A01P_DIFORTH==2 & A01P_RLGG==2 & A01P_DIFORA==2 & A01P_PBCOM==2 ~ "No",
                         .default = NA_character_)
    ) %>%
    
    
    #recode new variables
    mutate(
    # recode if stopped education in primary or secondary
    educ_2m=recode(educ_2m,
                   "0"="HighSchool_",
                   "1"="HighSchool_",
                   "2"="HighSchool_",
                   "3"="HighSchool_",
                   "4"="BacPlus2",
                   "5"="MoreThanBacPlus2",
                   
                     "-Inf"=NA_character_,
                     .default = as.character(educ_2m)),
    
    #recode M02_EXPTAB
    M02_EXPTAB=recode(M02_EXPTAB,
                      "2"="2_",
                      "3"="2_",
                      "4"="2_",
                      "5"="2_",
                      
                      "-Inf"=NA_character_,
                      .default = as.character(M02_EXPTAB)),
    
    #create a single tobacco variable from 
    #M00M2_TABAG = smoking during pregnancy
    #M00M2_EXPTAB = exposure to tobacco during pregnancy
    #M02_EXPTAB = exposure to tobacco at 2 months
    #
    # TOBACCO=case_when(M00M2_TABAG==1 & M02_EXPTAB=="1" ~ "Exp_PreNMum_NotExp_PostN",
    #                   M00M2_TABAG==0 & M00M2_EXPTAB=="Yes" & M02_EXPTAB=="1" ~ "Exp_PreNPass_NotExp_PostN",
    #                   M00M2_TABAG==1 & M02_EXPTAB=="2_" ~ "Exp_PreN_PostN",
    #                   M00M2_TABAG==0 & M00M2_EXPTAB=="Yes" & M02_EXPTAB=="2_" ~ "Exp_PreN_PostN",
    #                   M00M2_TABAG==0 & M00M2_EXPTAB=="No" & M02_EXPTAB=="2_" ~ "NotExp_PreN_Exp_PostN",
    #                   M00M2_TABAG==0 & M00M2_EXPTAB=="No" & M02_EXPTAB=="1" ~ "NotExp_PreN_PostN",
    #                   .default = NA_character_
    # ),
    
    TOBACCO=case_when(M00M2_TABAG==1  ~ "Exp_PreNMum",
                      M00M2_TABAG==0 & M00M2_EXPTAB=="Yes" ~ "Exp_PreNPass",
                      M00M2_TABAG==0 & M00M2_EXPTAB=="No" ~ "NotExp_PreN",
                      .default = NA_character_)
    )  %>% ungroup()   %>%
    
    #if you wanted to categorize BMI
    # M00M2_BMIMAVTG=case_when(M00M2_BMIMAVTG<=18.5 ~ "Low",
    #                          M00M2_BMIMAVTG>18.5 & M00M2_BMIMAVTG<=25 ~ "Normal",
    #                          M00M2_BMIMAVTG>25 & M00M2_BMIMAVTG<=30 ~ "High",
    #                          M00M2_BMIMAVTG>30 ~ "Obese"),
    #   
    # )
    
  
  #if you wanted to make family hx of learning difficulties a single variable grouping mum and dad 
  # A01_HxNDD=case_when(A01M_HxNDD=="Yes" | A01P_HxNDD=="Yes" ~ "Yes",
  #                     A01M_HxNDD=="No" & A01P_HxNDD=="No" ~ "No")
  
   mutate(
      M00M2_FQALCOOL=recode(M00M2_FQALCOOL, 
                            "0"="Never",
                            "1"="More",
                            "6"="More",
                            "More_1m"="More",
                            
                            .default=as.character(M00M2_FQALCOOL)),
      
      M00M3_FQCAFE=recode(M00M3_FQCAFE,
                          "6"="More_1_D",
                          .default="Less_1_D",
                          .missing = NA_character_),
      # # #0, <1/M, 1-3/M, 1/7, 2-5/7, 7/7, xxx/7
      #"M00M3_CAFE", #Never, More, Same, Less
      #"M00M3_FQCAFE", #6 categ
      
      M00M3_POISGEN=recode(M00M3_POISGEN,
                           "0"="Less_1m",
                           "1"="Less_1m",
                           "2"="More_1m",
                           "3"="More_1m",
                           "3"="More_1m",
                           "More_2_W"="More_1m"
                           ),
      #"M00M3_POISGEN",# Never, <1/31, 1-3/31, 1/7, 2-5/7, 7/7, xxx7/7, DNK
      # M00M3_CONSOPOIS More than once a month in the third trimester No, Yes
      
      M00M3_MGOMEGA3=recode(M00M3_MGOMEGA3, #passive other than home
                            "0"="No",
                            "1"="1",
                            "2"="2_",
                            "3"="2_",
                            "More_2_w"="2_",
                            
                            .default=NA_character_),
      #"M00M3_ACIDEG", # Yes
      #"M00M3_MGOMEGA3",# Never, <1/7, x/7, ~7/7, Always
      
      TOBACCO = recode(TOBACCO, 
                       "Exp_PreNMum"="Exp_PreN",
                       "Exp_PreNPass"="Exp_PreN",
                       "NotExp_PreN"="NotExp_PreN"),
      
      M00X4_TAU2010=recode(M00X4_TAU2010,
                           "0"="0_1",
                           "1"="0_1",
                           "2"="2_3_4",
                           "3"="2_3_4",
                           "4"="2_3_4",
                           "5"="5_6",
                           "6"="5_6"), 
      
      M02M_TYPALI=recode(M02M_TYPALI,#= breast, breast + bottle, bottle
                         "1"="breast",
                         "2"="bb",
                         "3"="bottle") 
      )


# Explore NA/distributions using summary/table/hist
str(merged_col_obs)
colSums(is.na(merged_col_obs))
map(names(merged_col_obs), function(var_name) {
  column <- merged_col_obs[[var_name]]
  if (!is.numeric(column)) {
    ggplot(data.frame(x = column), aes(x)) +
      geom_bar(binwidth = 1, fill = "dodgerblue", color = "white") +
      labs(title = var_name) +       # Add variable name as title
      theme_minimal()  } else {
    ggplot(data.frame(x = column), aes(x)) +
      geom_histogram(binwidth = 1, fill = "dodgerblue", color = "white") +
      labs(title = var_name) +       # Add variable name as title
      theme_minimal()
  }
})


# merged_gp ----
variable_set_2 <- c(
  
  #Outcome
  "ALL",
  
  #Exposure
  "preN_count_above_90", "preN_consec_count_above_90", "preN_count_above_95", "preN_consec_count_above_95",
  "postN1_count_above_90", "postN1_consec_count_above_90", "postN1_count_above_95", "postN1_consec_count_above_95",
  "postN2_count_above_90", "postN2_consec_count_above_90", "postN2_count_above_95", "postN2_consec_count_above_95",
  "postN3_count_above_90", "postN3_consec_count_above_90", "postN3_count_above_95", "postN3_consec_count_above_95",
  
  #Cold TDC
  "preN_count_below_10", "preN_consec_count_below_10", "preN_count_below_05", "preN_consec_count_below_05",
  "postN1_count_below_10", "postN1_consec_count_below_10", "postN1_count_below_05", "postN1_consec_count_below_05",
  "postN2_count_below_10", "postN2_consec_count_below_10", "postN2_count_below_05", "postN2_consec_count_below_05",
  "postN3_count_below_10", "postN3_consec_count_below_10", "postN3_count_below_05", "postN3_consec_count_below_05",
  
  # EDI, NDVI, Urban
  "EDI_2007","EDI_2011","EDI_2015",
  "buffer100m_ndvi_ete_2010","buffer100m_ndvi_ete_2011","buffer100m_ndvi_ete_2012","buffer100m_ndvi_ete_2013","buffer100m_ndvi_ete_2014", # 
  "M00X4_TAU2010",
  
  # CSP mum and dad
  "M00M2_CSP1M", 
  "M00M2_CSP1P",
  
  # From EQR12 (Socio-Demo)
  "c_emp_2m","c_emp_1y","c_emp_2y","c_emp_3y", #ceux qui cohabitent ont une activite pro 
  "sib_2m","sib_1y","sib_2y","sib_3y", #siblings
  "Child_hhld_2m","child_hhld_1y","child_hhld_2y","child_hhld_3y",#Ou vit l'enfant #"A02M_EFVIT",
  "revenu_part_qui_2m","revenu_part_qui_1y","revenu_part_qui_2y","revenu_part_qui_3y", # 
  "house_ownership_2m","house_ownership_1y","house_ownership_2y","house_ownership_3y", #
  
  # From EQR49 (Languages)
  "NBLANGMEN_2m", # "NBLANGMEN_1y",  "NBLANGMEN_2y",  "NBLANGMEN_3y", #nb of languages spoken at home
  
  # Moves
  "A01R_DEMENAG","A02R_DEMENAG","A03R_DEMENAG",
  
  #only one measure
  "educ_2m",  ### #higher of meduc & feduc,
  #"feduc_2m","meduc_2m",
  "imm", # migration status # 
  #"mimm","fimm",
  
  "M00M1_VAGUE", #season of birth
  #region
  
  "SEXE_ENF",
  
  "M00M2_AGEM", #mother's age at birth
  "M00M2_AGEP", #father's age at birth
  
  "M00M2_ENFGANT", 
  #"M00M2_FAF",
  
  "M00M2_BMIMAVTG",
  #"M00M2_POIMAVTG",
  #"M00M2_TAIM",
  
  "M00M2_FQALCOOL", #freq etOH (ordinal)
  
  "M00M3_FQCAFE", #6 categ
  
  "M00M3_POISGEN",# Never, <1/31, 1-3/31, 1/7, 2-5/7, 7/7, xxx7/7, always
  
   "M00M3_VITB9", # No Yes DNK
  #'M00M3_VITAB9',
  #'M00M2_ACIDEFOL',
  
  "M00M3_MGOMEGA3",# Never, <1/7, x/7, ~7/7, Allways
  
  #### 2M
  "M02M_TYPALI", # = breast, breast + bottle, bottle
  
  #tobacco
  "TOBACCO",
  #"M02_EXPTAB", #   
  #"M02M_EXPTAB", # exposure to tobacco 5 categ
  #"M02P_EXPTAB", #
  #"M00M2_TABAG",#  
  
  # "M00M2_EXPTAB", ##passive tobacco (home) 
  #"M00M2_EXPTABD", 
  #"M00M2_EXPTABLF", 
  
  #### 1Y
  "A01M_HxNDD",
  "A01P_HxNDD",
  #'A01M_DIFMATH',
  #'A01M_DIFLIR',
  #'A01M_DIFORTH',
  #'A01M_RLGG',
  #'A01M_DIFORA',
  #'A01M_PBCOM',
  #'A01P_DIFMATH',
  #'A01P_DIFLIR',
  #'A01P_DIFORTH',
  #'A01P_RLGG',
  #'A01P_DIFORA',
  #'A01P_PBCOM',
  
  #### 3Y
  "A03F_AGE3A"# age in months
  
)

#Group/arrange levels based on 30-2% per level & not too many levels & further steps
merged_gp <- merged_col_obs %>%
  select(all_of(variable_set_2),id_Dem806_585_GB)

char_N <- merged_gp %>%
  select(-id_Dem806_585_GB) %>%
  select(where(~ is.character(.))) %>%
  # Gather to long format to summarize levels
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  # Count levels including NA
  group_by(variable, value) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Calculate percentage within variable
  group_by(variable) %>%
  mutate(percentage = 100 * count / sum(count)) %>%
  ungroup() %>%
  # Replace NA level name for clarity
  mutate(level = ifelse(is.na(value), "NA", as.character(value))) %>%
  select(variable, level, count, percentage)


num_N <- merged_gp %>%
  select(-id_Dem806_585_GB) %>%
  select(where(~ (is.numeric(.) && n_distinct(.) < 10))) %>%
  # Gather to long format to summarize levels
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  # Count levels including NA
  group_by(variable, value) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Calculate percentage within variable
  group_by(variable) %>%
  mutate(percentage = 100 * count / sum(count)) %>%
  ungroup() %>%
  # Replace NA level name for clarity
  mutate(level = ifelse(is.na(value), "NA", as.character(value))) %>%
  select(variable, level, count, percentage)

# house_ownership_3y-Free
# c_emp_3y-None 138 1.5575621; # c_emp_2y-None 152 1.7155756
# educ_2m-_2 155 1.7494357



# plot var-outcome ----
# & biV, multV and redo merged_gp if necessary
# Exclude ID column and pivot data longer
cols_to_plot <- setdiff(names(merged_gp), c("id_Dem806_585_GB", "ALL"))

# Loop over each column to plot
p <- list()

for (col in cols_to_plot) {
  xvar <- merged_gp[[col]]
  
  if (is.character(xvar) || (is.numeric(xvar) && n_distinct(xvar) < 10)) {
    p[[col]] <- ggplot(merged_gp, aes(x = factor(.data[[col]]), y = ALL)) +
      geom_boxplot(alpha = 0.5) +
      #geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +
      labs(title = paste("ALL vs", col), x = col, y = "ALL") +
      theme_minimal()
  } else {
    p[[col]] <- ggplot(merged_gp, aes_string(x = col, y = "ALL")) +
      geom_point(alpha = 0.6) +
      geom_smooth(method=loess)+
      labs(title = paste("ALL vs", col), x = col, y = "ALL") +
      theme_minimal()
  }
}
p[["house_ownership_2m"]]; p[["house_ownership_1y"]];p[["house_ownership_2y"]]; p[["house_ownership_3y"]]#Free with Rent_HLM?
p[["c_emp_2m"]]; p[["c_emp_1y"]];p[["c_emp_2y"]]; p[["c_emp_3y"]]#prob do not concat?
p[["educ_2m"]]#can concat _2 and 3
#p

#merged_final ----
#Add on last set of variables (to make TDC)
merged_final <- merged_gp %>%
  mutate(M02R_DEMENAG="No") #%>%
  # mutate(EDI_2010=EDI_2007,
  #        EDI_2011=EDI_2011,
  #        EDI_2012= case_when(A02R_DEMENAG=="Yes" ~ NA_integer_,
  #                            TRUE ~ as.numeric(EDI_2011)
  #                            ),
  #        EDI_2013=case_when(A03R_DEMENAG=="Yes" ~ NA_integer_,
  #                           TRUE ~ as.numeric(EDI_2012)
  #        ),
  #        EDI_2014=EDI_2015) %>%
  # mutate(M00X4_TAU2011= case_when(A01R_DEMENAG=="Yes" ~ NA_character_,
  #                            TRUE ~ as.character(M00X4_TAU2010)
  #        ),
  #        M00X4_TAU2012= case_when(A02R_DEMENAG=="Yes" ~ NA_character_,
  #                                 TRUE ~ as.character(M00X4_TAU2011)
  #        ),
  #        M00X4_TAU2013= case_when(A03R_DEMENAG=="Yes" ~ NA_character_,
  #                                 TRUE ~ as.character(M00X4_TAU2012)
  #        ),
  #        M00X4_TAU2014= case_when(A03R_DEMENAG=="Yes" ~ NA_character_,
  #                                 TRUE ~ as.character(M00X4_TAU2013)
  #        )
  #        ) %>%
  # relocate(EDI_2011,.before = EDI_2012) %>%
  # relocate(EDI_2007,.before = EDI_2011) %>%
  # relocate(EDI_2015,.after = EDI_2014) %>%
  # 
  # relocate(M00X4_TAU2010,.before = M00X4_TAU2011) %>%
  # 
  # relocate(A01R_DEMENAG,.after = M00X4_TAU2014) %>%
  # relocate(A02R_DEMENAG,.after = A01R_DEMENAG) %>%
  # relocate(A03R_DEMENAG,.after = A02R_DEMENAG)


variable_final_set <- c(
  
  #Outcome
  "ALL",
  
  #Exposure
  "preN_count_above_90", "preN_consec_count_above_90", "preN_count_above_95", "preN_consec_count_above_95",
  "postN1_count_above_90", "postN1_consec_count_above_90", "postN1_count_above_95", "postN1_consec_count_above_95",
  "postN2_count_above_90", "postN2_consec_count_above_90", "postN2_count_above_95", "postN2_consec_count_above_95",
  "postN3_count_above_90", "postN3_consec_count_above_90", "postN3_count_above_95", "postN3_consec_count_above_95",
  
  #Cold TDC
  "preN_count_below_10", "preN_consec_count_below_10", "preN_count_below_05", "preN_consec_count_below_05",
  "postN1_count_below_10", "postN1_consec_count_below_10", "postN1_count_below_05", "postN1_consec_count_below_05",
  "postN2_count_below_10", "postN2_consec_count_below_10", "postN2_count_below_05", "postN2_consec_count_below_05",
  "postN3_count_below_10", "postN3_consec_count_below_10", "postN3_count_below_05", "postN3_consec_count_below_05",
  
  # EDI, NDVI, Urban
  "EDI_2007","EDI_2011","EDI_2015",
  "M00X4_TAU2010",
  "buffer100m_ndvi_ete_2010","buffer100m_ndvi_ete_2011","buffer100m_ndvi_ete_2012","buffer100m_ndvi_ete_2013","buffer100m_ndvi_ete_2014", # 
  
  # CSP mum and dad
  "M00M2_CSP1M", 
  "M00M2_CSP1P",
  
  # From EQR12 (Socio-Demo)
  "c_emp_2m","c_emp_1y","c_emp_2y","c_emp_3y", #ceux qui cohabitent ont une activite pro 
  "sib_2m","sib_1y","sib_2y","sib_3y", #siblings
  "Child_hhld_2m","child_hhld_1y","child_hhld_2y","child_hhld_3y",#Ou vit l'enfant #"A02M_EFVIT",
  "revenu_part_qui_2m","revenu_part_qui_1y","revenu_part_qui_2y","revenu_part_qui_3y", # 
  "house_ownership_2m","house_ownership_1y","house_ownership_2y","house_ownership_3y", #
  
  # Moves
  "M02R_DEMENAG","A01R_DEMENAG","A02R_DEMENAG","A03R_DEMENAG",
  
  # From EQR49 (Languages)
  "NBLANGMEN_2m", # "NBLANGMEN_1y",  "NBLANGMEN_2y",  "NBLANGMEN_3y", #nb of languages spoken at home
  
  #only one measure
  "educ_2m",  ### #higher of meduc & feduc,
  "imm", # migration status # 

  "M00M1_VAGUE", #season of birth
  #region
  
  "SEXE_ENF",
  
  "M00M2_AGEM", #mother's age at birth
  "M00M2_AGEP", #father's age at birth
  
  "M00M2_ENFGANT", 
  
  "M00M2_BMIMAVTG",
  
  "M00M2_FQALCOOL", #freq etOH (ordinal)
  
  "M00M3_FQCAFE", #6 categ
  
  "M00M3_POISGEN",# Never, <1/31, 1-3/31, 1/7, 2-5/7, 7/7, xxx7/7, always
  
  "M00M3_VITB9", # No Yes DNK
  
  "M00M3_MGOMEGA3",# Never, <1/7, x/7, ~7/7, Allways
  
  #### 2M
  "M02M_TYPALI", # = breast, breast + bottle, bottle
  
  #tobacco
  "TOBACCO",
  
  #### 1Y
  "A01M_HxNDD",
  "A01P_HxNDD",
  
  #### 3Y
  "A03F_AGE3A"# age in months
  
)

merged_final <- merged_final %>%
  select(all_of(variable_final_set), id_Dem806_585_GB)
merged_final$revenu_part_qui_2m <- as.character(merged_final$revenu_part_qui_2m)
merged_final$revenu_part_qui_1y <- as.character(merged_final$revenu_part_qui_1y)
merged_final$revenu_part_qui_2y <- as.character(merged_final$revenu_part_qui_2y)
merged_final$revenu_part_qui_3y <- as.character(merged_final$revenu_part_qui_3y)
merged_final$M00M1_VAGUE <- as.character(merged_final$M00M1_VAGUE)


#Add censoring ----
#censoring should occur after the last observed value
#then NA everywhere, except future C which should indicate 0
#Ask AI to code it for you

trt_suffix <- "consec_count_above_95"
trt <- c(paste0("preN_", trt_suffix),
         paste0("postN1_", trt_suffix),
         paste0("postN2_", trt_suffix),
         paste0("postN3_", trt_suffix)
)
cold_suffix <- "count_below_05"

baseline <- c(
  
  # EDI
  "EDI_2007",#"EDI_2011","EDI_2015",
  "M00X4_TAU2010",
  
  # CSP mum and dad
  "M00M2_CSP1M", 
  "M00M2_CSP1P",
  
  # From EQR49 (Languages)
  "NBLANGMEN_2m", # "NBLANGMEN_1y",  "NBLANGMEN_2y",  "NBLANGMEN_3y", #nb of languages spoken at home
  
  #only one measure
  "educ_2m",  ### #higher of meduc & feduc,
  "imm", # migration status # 
  
  "M00M1_VAGUE", #season of birth
  #region
  
  "SEXE_ENF",
  
  "M00M2_AGEM", #mother's age at birth
  "M00M2_AGEP", #father's age at birth
  
  "M00M2_ENFGANT", 
  
  "M00M2_BMIMAVTG",
  
  "M00M2_FQALCOOL", #freq etOH (ordinal)
  
  "M00M3_FQCAFE", #6 categ
  
  "M00M3_POISGEN",# Never, <1/31, 1-3/31, 1/7, 2-5/7, 7/7, xxx7/7, always
  
  "M00M3_VITB9", # No Yes DNK
  
  "M00M3_MGOMEGA3",# Never, <1/7, x/7, ~7/7, Allways
  
  #### 2M
  "M02M_TYPALI", # = breast, breast + bottle, bottle
  
  #tobacco
  "TOBACCO",
  
  #### 1Y
  "A01M_HxNDD",
  "A01P_HxNDD",
  
  #### 3Y
  "A03F_AGE3A"# age in months
)


time_vary <- list(
  
  preN=c(paste0("preN_",cold_suffix),"buffer100m_ndvi_ete_2010","c_emp_2m","sib_2m", "Child_hhld_2m","revenu_part_qui_2m","house_ownership_2m"),#"M02R_DEMENAG"
  postN1=c(paste0("postN1_",cold_suffix),"buffer100m_ndvi_ete_2011","c_emp_1y","sib_1y",  "child_hhld_1y","revenu_part_qui_1y","house_ownership_1y","A01R_DEMENAG"),
  postN2=c(paste0("postN2_",cold_suffix),"buffer100m_ndvi_ete_2012", "c_emp_2y","sib_2y",  "child_hhld_2y","revenu_part_qui_2y","house_ownership_2y","A02R_DEMENAG"),
  postN3=c(paste0("postN3_",cold_suffix),"buffer100m_ndvi_ete_2013","c_emp_3y","sib_3y",  "child_hhld_3y","revenu_part_qui_3y","house_ownership_3y","A03R_DEMENAG")
  
)


# 1. Define ids to be censored depending on time period
# baseline_vars <- NULL#baseline
preN_NotExp_id     <- sub("_.*", "",c(null_preN, all_names_NA_preN)) %>% unique %>% as.numeric# time_vary$preN)
postN1_NotExp_id   <- sub("_.*", "",c(null_postN, all_names_NA_postN)) %>% unique %>% as.numeric# time_vary$postN1)
postN2_NotExp_id   <- sub("_.*", "",c(null_postN, all_names_NA_postN)) %>% unique %>% as.numeric# time_vary$postN2)
postN3_NotExp_id   <- sub("_.*", "",c(null_postN, all_names_NA_postN)) %>% unique %>% as.numeric# time_vary$postN3)
# outcome_NotExp_id   <- "ALL"

# 2. Build the sequential, monotonic indicators
# C_0 is 1 only if BOTH baseline and (no multiple births, presence of preN_Exp)
merged_final$C_0 <- 1
merged_final$C_0[merged_final$id_Dem806_585_GB %in% id_MB$id_Dem806_585_GB] <- 0L
merged_final$C_0[merged_final$id_Dem806_585_GB %in% preN_NotExp_id] <- 0L
sum(merged_final$C_0)

# C_preN is 1 only if C_0 was 1 AND postN1 Exp is complete
merged_final$C_preN <- merged_final$C_0
merged_final$C_preN[merged_final$id_Dem806_585_GB %in% postN1_NotExp_id] <- 0L
sum(merged_final$C_preN)

# C_postN1 is 1 only if C_preN was 1 AND postN2 Exp are complete
merged_final$C_postN1 <- merged_final$C_preN

# C_postN2 is 1 only if C_postN1 was 1 AND postN3 Exp are complete
merged_final$C_postN2 <- merged_final$C_preN

# C_postN3 is 1 only if C_postN2 was 1 AND outcome variable is complete
merged_final$C_postN3 <- merged_final$C_preN

print(colSums(merged_final[, c("C_0", "C_preN", "C_postN1", "C_postN2", "C_postN3")]))


#merged_ignore ----
# 1. Build the ordered vector of column names
ordered_columns <- c(
  
  "id_Dem806_585_GB",
  
  # Step 0: Initial indicator
  "C_0",
  
  # Step 1: Baseline and preN block
  baseline,
  trt[1],
  time_vary$preN,
  "C_preN",

  # Step 2: postN1 block
  trt[2],
  time_vary$postN1,
  "C_postN1",
  
  # Step 3: postN2 block
  trt[3],
  time_vary$postN2,
  "C_postN2",
  
  # Step 4: postN3 block
  trt[4],
  time_vary$postN3,
  "C_postN3",
  
  #Outcome
  "ALL"
  
)

# 2. Create ordered dataframe
merged_ordered <- merged_final[, ordered_columns, drop = FALSE]

# 3. Impute all except missing Exp_preN and MB
merged_ignore.C_0 <- merged_ordered %>%
  filter(!id_Dem806_585_GB %in% c(pull(id_MB),preN_NotExp_id)) %>%
  select(id_Dem806_585_GB, all_of(c(baseline, trt[1], time_vary[[1]])))

merged_ignore_with_id_Cs <- merged_final %>%
  filter(!id_Dem806_585_GB %in% c(pull(id_MB),preN_NotExp_id)) 

merged_ignore <- merged_ignore_with_id_Cs %>%
  select(-id_Dem806_585_GB,-"C_0",-"C_preN",-"C_postN1",-"C_postN2",-"C_postN3")
  

# CHECK corr, naniar and drymice 
# RMV var/cases: obvious rmv (no value in observation) and more strategic rmv (influx-outflux)

#based on n_miss in individual cases and predefined threshold
#Function to better understand the distribution of missing values among variables and participants using naniar
try_missing <- function(data) {
  (miss_var_summary_obs <-naniar::miss_var_summary(data)); print(miss_var_summary_obs)
  (miss_var_table_obs <- naniar::miss_var_table(data) %>%
      mutate(pct_miss_in_var=n_miss_in_var*100/nrow(data), .after=n_miss_in_var))
  print(miss_var_table_obs)
  (miss_case_summary_obs <- naniar::miss_case_summary(data))
  print(miss_case_summary_obs)
  (miss_case_table_obs <- naniar::miss_case_table(data) %>%
      mutate(pct_miss_in_case=n_miss_in_case*100/ncol(data), .after=n_miss_in_case) %>%
      arrange(desc(pct_miss_in_case))
  )
  print(miss_case_table_obs)
  return(miss_case_summary_obs)
}

miss_case_summary_obs.C_0 <- try_missing(merged_ignore.C_0)

#Find a threshold criteria
p_miss_indiv <- 75
keep_id <- miss_case_summary_obs.C_0 %>%
  filter(pct_miss<=p_miss_indiv) %>%
  dplyr::select(case)
merged_ignore.C_0 <- merged_ignore.C_0[sort(keep_id$case),] 

md_pattern <- md.pattern(merged_ignore)
md.pairs(merged_ignore)
fx <- flux(merged_ignore);
plot(fx$influx, fx$outflux, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Influx", ylab = "Outflux", main = "Flux Plot")
text(fx$influx, fx$outflux, row.names(fx), pos = 4, cex = 0.8)

merged_ignore <- merged_ignore %>%
  mutate(across(
    .cols = where(function(x) is.character(x) || (is.numeric(x) && n_distinct(x) < 10)),
    .fns = as.factor
  ))
str(merged_ignore)
dryrun <- mice(merged_ignore, maxit = 0, print = FALSE)
# Inspect results of the dry run
print(dryrun$method)
print(dryrun$predictorMatrix)
print(dryrun$nmis)
print(dryrun$loggedEvents)
#Remove any constant and collinear variables before imputation.

# merged_imputed ----
# imp model, beware IA/non-linear, aux var, squeeze, post and passive imputation (trsf var e.g. BMI). sensitivity anal (MNAR). Data leak (ignore).
# remove id col
merged_imputed <- mice(merged_ignore, m=m, maxit = maxit, print = T)


# imputation dx ----
# (inc. Table1Imputed/NonImputed, density, strip) - warnings - logged events - FMI/LAMBDA ...
warnings()
merged_imputed$loggedEvents
colSums(is.na(merged_imputed %>% complete("long")))
colSums(is.na(merged_imputed %>% complete(1)))

plot(merged_imputed) #convergence
stripplot(merged_imputed)#values for imputed datasets imputed and non-imputed points
densityplot(merged_imputed)


#savings ----
#add id to merged_final
merged_imputed$id <- merged_ignore_with_id_Cs$id_Dem806_585_GB
save(merged_ignore_with_id_Cs, merged_imputed, file = paste0(Exp, "_merged_imputed_IPCW.RData" ))
