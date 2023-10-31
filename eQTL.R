
#???ð?
library(TwoSampleMR)


#exposureID="ebi-a-GCST004603"         #??¶????ID
outcomeID="ukb-a-87"   
#load("FBL.rdata")

#??ȡ??¶????
#exposure_dat <- exp_dat
exposure_dat <- read.table("gene.txt", sep="\t", header=TRUE)
# 完善数据
exposure_dat$eaf.exposure <- NA
exposure_dat$samplesize.exposure <- NA
exposure_dat$pos.exposure <- exposure_dat$BP 
exposure_dat$mr_keep.exposure <- NA
exposure_dat$pval_origin.exposure <- NA
exposure_dat$f <- NA
library(dplyr)

exposure_dat <- exposure_dat %>%
  rename(
    id.exposure = SNP,
    effect_allele.exposure = A1,
    other_allele.exposure = A2,
    beta.exposure = b,
    se.exposure = SE
  )
exposure_dat <- exposure_dat %>%
  rename(
    SNP = id.exposure
  ) %>%
  mutate(
    exposure = Gene
  )

exposure_dat$id.exposure <- exposure_dat$exposure


#??ȡ????????
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)

head(exposure_dat)
head(outcome_dat)
#????¶???ݺͽ??????ݺϲ?
dat <- harmonise_data(exposure_dat, outcome_dat)

#?????????ϵ¶????????Ĺ??߱?��
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

#?ϵ¶???????????
mrResult=mr(dat)
#ѡ???ϵ¶????????ķ???
#mr_method_list()$obj
#mrResult=mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
#?Խ???????ORֵ?ļ???
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#?????Է???
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#??Ч?Լ???
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#????ɢ??ͼ
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#ɭ??ͼ
res_single=mr_singlesnp(dat)      #?õ?ÿ?????߱?��?Խ??ֵ?Ӱ??
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#©??ͼ
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#??һ???????Է???
pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()


######Vsave.image("WORK.Rdata")