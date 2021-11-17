setwd("/Users/namuhanz/cancer/transcripts/enrichment analysis/file")
library(veccompare)

bf_df <- read.csv("./bladder_female_msigdb.csv", header = T)
bm_df <- read.csv("./bladder_male_msigdb.csv", header = T)
cf_df <- read.csv("./colon_female_msigdb.csv", header = T)
cm_df <- read.csv("./colon_male_msigdb.csv", header = T)
kf_df <- read.csv("./kidney_female_msigdb.csv", header = T)
km_df <- read.csv("./kidney_male_msigdb.csv", header = T)
lf_df <- read.csv("./liver_female_msigdb.csv", header = T)
lm_df <- read.csv("./liver_male_msigdb.csv", header = T)
sf_df <- read.csv("./stomach_female_msigdb.csv", header = T)
sm_df <- read.csv("./stomach_male_msigdb.csv", header = T)
tf_df <- read.csv("./thyroid_female_msigdb.csv", header = T)
tm_df <- read.csv("./thyroid_male_msigdb.csv", header = T)

bf <- bf_df[bf_df$p.adjust<=0.05,]$ID
bm <- bm_df[bm_df$p.adjust<=0.05,]$ID
cf <- cf_df[cf_df$p.adjust<=0.05,]$ID
cm <- cm_df[cm_df$p.adjust<=0.05,]$ID
kf <- kf_df[kf_df$p.adjust<=0.05,]$ID
km <- km_df[km_df$p.adjust<=0.05,]$ID
lf <- lf_df[lf_df$p.adjust<=0.05,]$ID
lm <- lm_df[lm_df$p.adjust<=0.05,]$ID
sf <- sf_df[sf_df$p.adjust<=0.05,]$ID
sm <- sm_df[sm_df$p.adjust<=0.05,]$ID
tf <- tf_df[tf_df$p.adjust<=0.05,]$ID
tm <- tm_df[tm_df$p.adjust<=0.05,]$ID

genesetlist <- list(bf,bm,cf,cm,kf,km,lf,lm,sf,sm,tf,tm)
genesetlist_female <- list(bf,cf,kf,lf,sf,tf)
genesetlist_male <- list(bm,cm,km,lm,sm,tm)

universal_geneset <- Reduce(intersect, genesetlist)
universal_geneset_female <- Reduce(intersect, genesetlist_female)
universal_geneset_male <- Reduce(intersect, genesetlist_male)

write.csv(universal_geneset, file = "./universal gene sets.csv")
write.csv(universal_geneset_female, file = "./universal gene sets female.csv")
write.csv(universal_geneset_male, file = "./universal gene sets male.csv")

