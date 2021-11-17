setwd("/Users/namuhanz/cancer/transcripts/compare DE/file")
df <- read.table("./c2.xls", header = T)
dm <- read.table("./c3.xls", header = T)
fgl <- read.csv("./female_specific.csv", header = T)
mgl <- read.csv("./male_specific.csv", header = T)
nrf <- nrow(fgl)
nrm <- nrow(mgl)

female_genelist <- data.frame(gene_id = character(nrf),
                              log2FC = numeric(nrf))
male_genelist <- data.frame(gene_id = character(nrm),
                            log2FC = numeric(nrm))

female_genelist[,1] <- fgl[,1]
male_genelist[,1] <- mgl[,1]

for(gene in female_genelist$gene_id){
  female_genelist[female_genelist$gene_id==gene,2] <- df[df$gene_id==gene, "log2FoldChange"]
}

for(i in male_genelist$gene_id){
  male_genelist[male_genelist$gene_id==i,2] <- dm[dm$gene_id==i, "log2FoldChange"]
}

write.csv(female_genelist, file = "./female_genelist.csv")
write.csv(male_genelist, file = "./male_genelist.csv")

