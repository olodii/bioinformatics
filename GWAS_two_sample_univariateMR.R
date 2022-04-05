#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

#----------------------------------------------------------------------------------------------------------------#
# univariate MR: exposure와 health outcome 과의 연관성
#----------------------------------------------------------------------------------------------------------------#
setwd("C://Users//USER//Desktop//bioinformatics")

library(TwoSampleMR)
library(dplyr)
library(ggplot2)
select <- dplyr::select

ao <- available_outcomes()

### exposure
exp_id <- c('SNP_id') 
exposure_dat <- extract_instruments(exp_id)

### outcome
oc_id <- c('outcome1_from_consortium',
           'outcome2_from_consortium',
           'outcome3_from_consortium6')

outcome_dat <- extract_outcome_data(exposure_dat$SNP, oc_id)

### clumping
exposure_dat <- clump_data(exposure_dat)

### analysis
dat <- harmonise_data(exposure_dat, outcome_dat)

### results
res <- mr(dat)
res %>% View()

### OR, 95%CI 포함
odds_res <- generate_odds_ratios(res)

### representative method: focus on IVW method
sub_odds_res <- subset(odds_res, method=="Inverse variance weighted") %>% select(method,outcome,nsnp,or,or_lci95,or_uci95,pval)


names(sub_odds_res) <- c("method","outcome","nsnp","OR","lCI_95","uCI_95","pval")

sub_odds_res$OR <- round(sub_odds_res$OR, 3)
sub_odds_res$lCI_95 <- round(sub_odds_res$lCI_95, 3)
sub_odds_res$uCI_95 <- round(sub_odds_res$uCI_95, 3)
sub_odds_res$pval <- round(sub_odds_res$pval, 3)

res
res$sig <- ifelse(res$pval<0.001,"***",ifelse(res$pval<0.01,"**",ifelse(res$pval<0.05,"*","")))

names(res)
res %>% select(exposure,outcome,method, b,se,pval,sig,nsnp) %>% View()

### Sensitivity: Horizontal pleiotropy examined through MR-Egger regression
mr_pleiotropy_test(dat) %>% View()
pleio <- mr_pleiotropy_test(dat)
pleio$egger_intercept.r <- round(pleio$egger_intercept,digits=3)

### scatter plot : Analysing 'SNP_id' on 'outcomes_from_consortium'
p1 <- mr_scatter_plot(res, dat)
x11();p1[[1]]
x11();p1[[2]]
x11();p1[[3]]

# if you want to save image
#ggsave(p1[[1]], file="scatterplot_outcome1.png", width=7, height=7)
#ggsave(p1[[2]], file="scatterplot_outcome2.png", width=7, height=7)

### forest plot
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
x11();p2[[1]]
x11();p2[[2]]
x11();p2[[3]]

# if you want to save image
#ggsave(p2[[1]], file="forestplot_outcome1.png", width=7, height=7)
#ggsave(p2[[2]], file="forestplot_outcome2.png", width=7, height=7)

### leave-one-out plot
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
x11();p3[[1]]
x11();p3[[2]]
x11();p3[[3]]

# if you want to save image
#ggsave(p3[[1]], file="leavoneoutplot_outcome1.png", width=7, height=7)
#ggsave(p3[[2]], file="leavoneoutplot_outcome2.png", width=7, height=7)

### funnel plot
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
x11();p4[[1]]

# if you want to save image
#ggsave(p4[[1]], file="funnelplot_outcome1.png", width=7, height=7)
#ggsave(p4[[2]], file="funnelplot_outcome2.png", width=7, height=7)

### REPORT
mr_report(dat)


