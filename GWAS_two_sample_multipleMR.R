
#----------------------------------------------------------------------------------------------------------------#
# Multiple MR : exposure와 여러 outcomes/여러 exposure와 outcome1간의 연관성 확인
#----------------------------------------------------------------------------------------------------------------#
library(TwoSampleMR)
library(dplyr)
library(ggplot2)
select <- dplyr::select

ao <- available_outcomes()

### exposure and outcome
exp_id <- c('SNP_id') # main exposure


exposure_dat <- extract_instruments(exp_id)

oc_id <- c('ukb-b-10753' ,# Diabetes - ex.
           'ukb-b-8714') # Stroke - ex.

outcome_dat <- extract_outcome_data(exposure_dat$SNP, oc_id)

### clumping
exposure_dat <- clump_data(exposure_dat)

### analysis
dat <- harmonise_data(exposure_dat, outcome_dat)

### results
res <- mr(dat)

res$sig <- ifelse(res$pval<0.001,"***",ifelse(res$pval<0.01,"**",ifelse(res$pval<0.05,"*","")))
res %>% select(exposure,outcome,method, b,se,pval,sig,nsnp) %>% View()


