setwd("/Users/namuhanz/cancer/mRNA-methylation/file")

femaledm_hyper_target <- read.csv("./femaledm_hyper_target.csv", header = T)
femaledm_hypo_target <- read.csv("./femaledm_hypo_target.csv", header = T)
maledm_hyper_target <- read.csv("./maledm_hyper_target.csv", header = T)
maledm_hypo_target <- read.csv("./maledm_hypo_target.csv", header = T)

female_hyper_target <- femaledm_hyper_target$x
female_hypo_target <- femaledm_hypo_target$x
male_hyper_target <- maledm_hyper_target$x
male_hypo_target <- maledm_hypo_target$x

femalede <- read.csv("./female_specific.csv", header = T)
malede <- read.csv("./male_specific.csv", header = T)

female_up <- femalede[femalede$sig=="Up",]$gene_id
female_down <- femalede[femalede$sig=="Down",]$gene_id
male_up <- malede[malede$sig=="Up",]$gene_id
male_down <- malede[malede$sig=="Down",]$gene_id

femaletrueup <- intersect(female_up,female_hypo_target)
femaletruedown <- intersect(female_down,female_hyper_target)
maletrueup <- intersect(male_up,male_hypo_target)
maletruedown <- intersect(male_down,male_hyper_target)

write.csv(femaletrueup, file = "./femaledm_trueup.csv")
write.csv(femaletruedown, file = "./femaledm_truedown.csv")
write.csv(maletrueup, file = "./maledm_trueup.csv")
write.csv(maletruedown, file = "./maledm_truedown.csv")


print(c(length(femaletrueup),length(femaletruedown),length(maletrueup),length(maletruedown)))

























