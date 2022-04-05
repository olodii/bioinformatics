
#-----------------------HoW to Post Results?----------------------------------------------------#


head(res)

res<-mr(dat)

split_outcome(res)

generate_odds_ratios(res)

subset_on_method(res)



res<-mr(dat)
het<-mr_heterogeneity(dat)
plt<-mr_pleiotropy_test(dat)
sin<-mr_singlesnp(dat)
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,Exp=T,split.exposure=F,split.outcome=T)
head(all_res[,c("Method","outcome","exposure","nsnp","b","se","pval","intercept","intercept_se","intercept_pval","Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])