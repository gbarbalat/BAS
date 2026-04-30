rm(list=ls())

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mice)
library(data.table)
library(stringr)
library(lmtp)
library(glmnet)

data_path <- "/bettik/barbalag/BAS/"

debug <- FALSE
Exp <- "Tmin"#Tmean, Tmax, Tmin
trt_suffix <- "consec_count_above_95" #all_trt_suffix <- c("count_above_90", "count_above_95", "consec_count_above_95", "consec_count_above_90") ; 
use_Exp_pct <- TRUE # use pct of heat per period or raw number of heat event per period
use_shift <- FALSE # use shift vs. shifted argument in lmtp functions
m <- 10 #imputed dataset
folds <- 5 #according to Philipps et al. 2022, for this N, needs to be between 5 and 10.

all_trim <- c(0.90, 0.925, 0.95, 0.975); #c(0.90, 0.925, 0.95, 0.975, 0.99, 1); 
all_CF <- c("CF1", "CF2", "CF3") 
all_period <- c("preN", "postN1", "postN2", "postN3")

if (debug) {
learners <- c("SL.mean")
args <- c(1,1,1)
} else {
learners <- c("SL.mean","SL.glm", "SL.glmnet", "SL.earth", "SL.ranger", "SL.xgboost")
args <- commandArgs(TRUE); 
}
trim <- all_trim[as.numeric(args[1])]; print(trim)
CF <- all_CF[as.numeric(args[2])]; print(trt_suffix)
period <- all_period[as.numeric(args[3])]; print(period)

if (trt_suffix == "count_above_95") {
  cold_suffix <- "count_below_05"
} else if (trt_suffix == "consec_count_above_95") {
  cold_suffix <- "consec_count_below_05"
} else if (trt_suffix == "count_above_90") {
  cold_suffix <- "count_below_10"
} else if (trt_suffix == "consec_count_above_90") {
  cold_suffix <- "consec_count_below_10"
}


trt <- c(paste0("preN_", trt_suffix),
         paste0("postN1_", trt_suffix),
         paste0("postN2_", trt_suffix),
         paste0("postN3_", trt_suffix)
         )


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


moderators <- c(
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
  
  "SEXE_ENF",
  
  "M00M2_AGEM", #mother's age at birth
  "M00M2_AGEP", #father's age at birth
  
  "M00M2_ENFGANT", 
  
  "M00M2_BMIMAVTG",
  
  #### 1Y
  "A01M_HxNDD",
  "A01P_HxNDD"
)



time_vary <- list(
 
  preN=c(paste0("preN_",cold_suffix),"buffer100m_ndvi_ete_2010","c_emp_2m","sib_2m", "Child_hhld_2m","revenu_part_qui_2m","house_ownership_2m"),#"M02R_DEMENAG"
  postN1=c(paste0("postN1_",cold_suffix),"buffer100m_ndvi_ete_2011","c_emp_1y","sib_1y",  "child_hhld_1y","revenu_part_qui_1y","house_ownership_1y","A01R_DEMENAG"),
  postN2=c(paste0("postN2_",cold_suffix),"buffer100m_ndvi_ete_2012", "c_emp_2y","sib_2y",  "child_hhld_2y","revenu_part_qui_2y","house_ownership_2y","A02R_DEMENAG"),
  postN3=c(paste0("postN3_",cold_suffix),"buffer100m_ndvi_ete_2013","c_emp_3y","sib_3y",  "child_hhld_3y","revenu_part_qui_3y","house_ownership_3y","A03R_DEMENAG")

)

load (file = paste0(data_path ,Exp, "_merged_imputed.RData" ))
load(file=paste0(data_path ,Exp, "_Exp_", CF, ".RData"))
load(file=paste0(data_path ,Exp, "_length_Exp.RData")) #length_Exp_preN #length_Exp_postN3

#if we want to use pct of days with Exposure higher than a threshold rather than number of days
names(length_Exp_preN) <- sub("_.*", "", names(length_Exp_preN))
length_Exp_preN <- length_Exp_preN[!duplicated(names(length_Exp_preN))]

names(length_Exp_postN3) <- sub("_.*", "", names(length_Exp_postN3))
length_Exp_postN3 <- length_Exp_postN3[!duplicated(names(length_Exp_postN3))]

Exp_as_pct <- function(which_df) {
  which_df <- which_df %>%
    mutate(
      # divide preN_ columns by length_Exp_preN[id]
      across(starts_with("preN_"),
             ~ .x / length_Exp_preN[as.character(id)]),
      
      # divide postN1_ columns by 365
      across(starts_with("postN1_"),
             ~ .x / 365),
      
      # divide postN2_ columns by 366
      across(starts_with("postN2_"),
             ~ .x / 366),
      
      # divide postN3_ columns by (length_Exp_postN3[id])
      across(starts_with("postN3_"),
             ~ .x / (length_Exp_postN3[as.character(id)]))
    )
}


#fn to run lmtp
run_lmtp_sdr <- function(data, mtp, shift, shifted) {

lmtp_sdr(
  data=data,
  trt=trt,
  outcome="ALL",
  baseline = baseline,
  time_vary = time_vary, #time_vary,#NULL
  cens = NULL,
  shift = shift,
  shifted = shifted,
  k = Inf,
  mtp = mtp,
  outcome_type = "continuous",
  id = NULL,
  bounds = NULL,
  learners_outcome = learners,#learners,"SL.glm"
  learners_trt = learners,#learners,"SL.glm"
  folds = folds,
  weights = NULL,
  control = lmtp_control(.trim=trim,
                         .learners_outcome_folds = folds,
                         .learners_trt_folds = folds)
)
}


#main lmtp fn to be run for each imputed dataset
run_lmtp_imputed_df <- function(i) {

  #factual df
  merged_imputed_thisone <- merged_imputed %>% complete(i)
  merged_imputed_thisone$id <- merged_imputed$id
  
  #CF df
  merged_imputed_CF <- merged_imputed_thisone;   merged_imputed_CF[paste0(period, "_", trt_suffix)] <- get(paste0("Exp_",CF))[paste0(period, "_", trt_suffix)]
  #merged_imputed_CF2 <- merged_imputed_thisone;   merged_imputed_CF2[paste0(period, "_", trt_suffix)] <- Exp_CF2[paste0(period, "_", trt_suffix)]  
  
if (use_Exp_pct) {
  merged_imputed_thisone <- Exp_as_pct(merged_imputed_thisone) 
  merged_imputed_CF <- Exp_as_pct(merged_imputed_CF) 
}

if (use_shift) {

if (CF=="CF1") {
	policy <- function(data, x) {
   	data[[x]] + 0.05
	}
} else if (CF=="CF2") {
	policy <- function(data, x) {
   	data[[x]] + 0.1
	}
} else if (CF=="CF3") {
	policy <- function(data, x) {
   	data[[x]] + 0.15
	}
}

set.seed(123); factual <- run_lmtp_sdr(data=merged_imputed_thisone,mtp= FALSE, shift=NULL, shifted=NULL)
set.seed(123); cfactual <- run_lmtp_sdr(data=merged_imputed_thisone,mtp=TRUE, shift=policy, shifted=NULL)
} else {

set.seed(123); factual <- run_lmtp_sdr(data=merged_imputed_thisone,mtp= FALSE, shift=NULL, shifted=NULL)
set.seed(123); cfactual <- run_lmtp_sdr(data=merged_imputed_thisone,mtp=TRUE, shift=NULL, shifted=merged_imputed_CF)
}
  

#standard lmtp contrast
contrast <- lmtp_contrast(cfactual, ref=factual)

print(i)

return(list(contrast=contrast,
            cfactual=cfactual,
            factual=factual, 
merged_imputed_thisone=merged_imputed_thisone,
merged_imputed_CF=merged_imputed_CF
            )
)

}


results_lmtp <- lapply(1:m, run_lmtp_imputed_df)
save(results_lmtp, file=paste0(Exp, "_", CF , "_trim", trim, "_", paste0(period,collapse="_"), "_", trt_suffix, ".RData"))