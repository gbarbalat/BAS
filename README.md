# BAS

## Paper  
submitted  

## Scripts  
- Pre-process_data_BAS.R  
- Process_data_BAS.R
- Process_data_BAS_CF.R does the same but for counterfactual exposures
- compare_inc_full.R and compare_imp_non_imp.R produce tables comparing full cohort and selected cohort as well as final sample imputed vs. not imputed. Uses zombie_process_data_for_full_sample.R as a helper.  
- anal_BAS_lmtp_CF.R runs the analysis
- anal_BAS_lmtp_CF_byPeriods.R runs the analysis for each period of analysis
- anal_results_DR.R analyses the density ratios
- anal_results_main.R performs the main plots
- anal_results_EffMdif.R performs the effect modification analysis


## Asked by reviewers
- Exclude premature babies and add Air Pollution as a co-exposure: use anal_BAS_lmtp_prema_Poll.R and Process_data_BAS_CoExp.R  
- Use alternative thresholds for Heat (already done but not in the paper)
- 3 months-windows during pregnancy: use Process_data_3monthsPregnancy_BAS.R which calls create_db_policy_3monthsPregnancy located in Process_data_3monthsPregnancy_BAS_CF.R. Also uses anal_BAS_lmtp_CF_byPeriods_3monthsPregnancy.R
- Use IPCW: use Process_data_BAS_IPCW.R then anal_BAS_lmtp_CF_IPCW.R


## sh scripts  
Run multiple scripts with each one having a full node  
In shell:  oarsub –S ./RUN_anal_lmtp_Prema_Poll.sh --array-param-file RUN_params.txt  
Job is launched as many times as they are rows in RUN_params.txt  
In .sh file (assuming two columns), last row: Rscript %name of script%  $1 $2  
In .R script:  
args=commandArgs(trailingOnly=TRUE)  
i=as.numeric(args[1]); j=as.numeric(args[2]);  
In RUN_paramts.txt your values of i and j  
1 2  
2 4  
3 6  

- RUN_anal_lmtp_CF_Prema_Poll.sh with anal_BAS_lmtp_CF_Prema_Poll.R $1 $2 and trim_CF.txt
- RUN_anal_lmtp_CF_byPeriods_3monthsPregnancyl.sh with anal_BAS_lmtp_CF_byPeriods_3monthsPregnancy.R $1 $2 $3 and trim_CF1_period.txt, trim_CF2_period.txt, trim_CF3_period.txt
- RUN_Process_BAS_IPCW.sh with Process_data_BAS_IPCW.R
- RUN_anal_lmtp_CF_IPCW.sh with anal_BAS_lmtp_CF_IPCW.R $1 $2 and trim_CF.txt
