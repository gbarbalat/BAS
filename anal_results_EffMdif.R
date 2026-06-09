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
  grepl("consec_count_above_95", rdata_files) &
  !grepl("trim1", rdata_files) &
  !grepl("trim0.99", rdata_files)
]

#for lmtp
out <- vector("list", length(rdata_files))
summary_df <- NULL

#for effect modification
pooled_results <- list()
glmnet_coeff_list <- list()# to store glmnet coeff by imputed datasets
Exp <- "Tmin"
data_path <- "/bettik/barbalag/BAS/"

for (f in rdata_files) {
  
  load(f)
  
#load("Tmin_CF3_trim0.9_trt_suffix_count_above_95.RData")
# results_lmtp[[1]]$contrast$vals$theta
# results_lmtp[[1]]$contrast$vals$std.error

# results_lmtp[[1]]$contrast$eifs %>% head
# (results_lmtp[[1]]$factual$eif - results_lmtp[[1]]$cfactual$eif) %>% head


# part I - pooled effect estimates see anal_results_Main.R ----

# partII - effect modification see anal_results_effect_mdif.R ----
#regress contrasted EIF onto baseline covariates
#23 baseline
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

#21 moderators
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
  
  #"M00M1_VAGUE", #season of birth
  
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
  #"M02M_TYPALI", # = breast, breast + bottle, bottle
  
  #tobacco
  "TOBACCO",
  

  #### 1Y
  "A01M_HxNDD",
  "A01P_HxNDD"
)
trt_suffix <- "consec_count_above_95" #all_trt_suffix <- c("count_above_90", "count_above_95", "consec_count_above_95", "consec_count_above_90") ; 
if (trt_suffix == "count_above_95") {
  cold_suffix <- "count_below_05"
} else if (trt_suffix == "consec_count_above_95") {
  cold_suffix <- "consec_count_below_05"
} else if (trt_suffix == "count_above_90") {
  cold_suffix <- "count_below_10"
} else if (trt_suffix == "consec_count_above_90") {
  cold_suffix <- "consec_count_below_10"
}
time_vary <- list(
 
  preN=c(paste0("preN_",cold_suffix),"buffer100m_ndvi_ete_2010","c_emp_2m","sib_2m", "Child_hhld_2m","revenu_part_qui_2m","house_ownership_2m"),#"M02R_DEMENAG"
  postN1=c(paste0("postN1_",cold_suffix),"buffer100m_ndvi_ete_2011","c_emp_1y","sib_1y",  "child_hhld_1y","revenu_part_qui_1y","house_ownership_1y","A01R_DEMENAG"),
  postN2=c(paste0("postN2_",cold_suffix),"buffer100m_ndvi_ete_2012", "c_emp_2y","sib_2y",  "child_hhld_2y","revenu_part_qui_2y","house_ownership_2y","A02R_DEMENAG"),
  postN3=c(paste0("postN3_",cold_suffix),"buffer100m_ndvi_ete_2013","c_emp_3y","sib_3y",  "child_hhld_3y","revenu_part_qui_3y","house_ownership_3y","A03R_DEMENAG")

)


# fn for effect modification on eif
fn_effect_mdif <- function(eif, merged_imputed_thisone) {

  merged_imputed_thisone$M02M_TYPALI <- relevel(factor(merged_imputed_thisone$M02M_TYPALI), ref = "bottle")  

  #glm
  predictors <- paste(merged_imputed_thisone %>% dplyr::select(all_of(moderators), unlist(time_vary$preN) %>% as.vector) %>% colnames, collapse = "+")
  glm_fit <- glm(formula=formula(paste("eif ~ ", predictors)), data=merged_imputed_thisone)
  theta <- coef(glm_fit)
  se <- summary(glm_fit)$coefficients[, "Std. Error"]
  
  
  #glmnet
  mm <- model.matrix(~ . -1, merged_imputed_thisone %>% dplyr::select(all_of(moderators), unlist(time_vary$preN) %>% as.vector))
  cv.glmnet_model <- cv.glmnet(y=eif, 
                               x=mm, 
                               alpha=1, #lasso(alpha=1) penalty 
                               nfolds=20)
  glmnet_fit <- glmnet(x=mm,
                       y=eif,  
                       family = "gaussian",  
                       alpha=1, #lasso penalty 
                       lambda = cv.glmnet_model$lambda.min )
  r2_range <- glmnet_fit$dev.ratio;
  glmnet_coeff <- as.matrix(coef(glmnet_fit,s=cv.glmnet_model$lambda.min ))
  
  #return(list(glm_fit=glm_fit, glmnet_fit=glmnet_fit)
  return(list(theta=theta, std.error=se, glmnet_coeff=glmnet_coeff))
  
}
load (file = paste0(data_path, Exp, "_merged_imputed.RData" ))

#loop over imputations
M <- length(results_lmtp) -> m # number of imputations

# Nb predictors
coef_names <- model.matrix(~., merged_imputed %>% complete(1) %>% dplyr::select(all_of(moderators), unlist(time_vary$preN) %>% as.vector)) %>% colnames

p <- length(coef_names)

#initialize theta and se matrices
theta_mat <- matrix(NA, nrow = p, ncol = m,
                    dimnames = list(coef_names, paste0("imp", 1:m)))

se_mat <- matrix(NA, nrow = p, ncol = m,
                 dimnames = list(coef_names, paste0("imp", 1:m)))

glmnet_coeff_mat <- NULL

for (i in seq_len(M)) {
  # Extract completed imputation
  merged_imputed_thisone <- merged_imputed %>% complete(i)
  
  # Compute your effect
  effect_mdif_i <- fn_effect_mdif(
    results_lmtp[[i]]$contrast$eifs,
    merged_imputed_thisone
  )
  
  theta_mat[, i] <- effect_mdif_i$theta
  se_mat[, i]    <- effect_mdif_i$std.error
  glmnet_coeff_mat    <- cbind(glmnet_coeff_mat , effect_mdif_i$glmnet_coeff)

  print(i)
}

#Rubin rule with df
rubin_pool_matrix <- function(theta_mat, se_mat) {
  
  stopifnot(is.matrix(theta_mat),
            is.matrix(se_mat),
            dim(theta_mat) == dim(se_mat))
  
  M <- ncol(theta_mat)
  p <- nrow(theta_mat)
  
  pooled_list <- lapply(seq_len(p), function(j) {
    
    theta_j <- theta_mat[j, ]
    se_j    <- se_mat[j, ]
    
    theta_bar <- mean(theta_j)
    U_bar <- mean(se_j^2)
    B <- var(theta_j)
    
    T_var <- U_bar + (1 + 1/M) * B
    se_pooled <- sqrt(T_var)
    
    # Barnard–Rubin degrees of freedom
    lambda <- ((1 + 1/M) * B) / T_var
    df <- (M - 1) / lambda^2
    
    t_stat <- theta_bar / se_pooled
    p_val <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
    
    ci_low  <- theta_bar - qt(0.975, df) * se_pooled
    ci_high <- theta_bar + qt(0.975, df) * se_pooled
    
    data.frame(
      estimate = theta_bar,
      se = se_pooled,
      df = df,
      p_value = p_val,
      ci_low = ci_low,
      ci_high = ci_high
    )
  })
  
  pooled_df <- do.call(rbind, pooled_list)
  pooled_df$term <- rownames(theta_mat)
  
  pooled_df %>% relocate(term)
}
pooled_results[[basename(f)]] <- rubin_pool_matrix(theta_mat, se_mat)
glmnet_coeff_list[[basename(f)]] <- glmnet_coeff_mat

print(basename(f))

}

pooled_results
save(glmnet_coeff_list,pooled_results, file="BAS_lmtp_EffMdif.RData")