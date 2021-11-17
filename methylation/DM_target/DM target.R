setwd("/Users/namuhanz/cancer/methylation/DM target/file")
dm <- read.csv("DMR_results_sample.type_Primary.Tumor_Solid.Tissue.Normal_pcut_0.05_meancut_0.25.csv", header = T)
dm_hyper <- dm[dm$status=="Hypermethylated in Primary Tumor",]
dm_hypo <- dm[dm$status=="Hypomethylated in Primary Tumor",]
hyper_cg <- dm_hyper$X
hypo_cg <- dm_hypo$X
cg_convert <- read.csv("../ref.csv", header = T)
rownames(cg_convert) <- cg_convert[,"cg_id"]
cg_convert <- cg_convert[,-1]

hyper_target <- character()
for (i in 1:length(hyper_cg)) {
  cg <- hyper_cg[i]
  hyper_target <- c(hyper_target,cg_convert[cg,]$gene_symbol)
}
hyper_target <- hyper_target[hyper_target!="."]
hyper_target <- unlist(strsplit(hyper_target, split=", "))
hyper_target <- unique(hyper_target)

hypo_target <- character()
for (i in 1:length(hypo_cg)) {
  cg <- hypo_cg[i]
  hypo_target <- c(hypo_target,cg_convert[cg,]$gene_symbol)
}

hypo_target <- hypo_target[hypo_target!="."]
hypo_target <- unlist(strsplit(hypo_target, split=", "))
hypo_target <- unique(hypo_target)

write.csv(hyper_target, file = "hyper-methyl targets.csv")
write.csv(hypo_target, file = "hypo-methyl targets.csv")

write.csv(hyper_cg, file = "hyper-methyl cg.csv")
write.csv(hypo_cg, file = "hypo-methyl cg.csv")










