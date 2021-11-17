setwd("/Users/namuhanz/cancer/mRNA-methylation/file")
# library(TCGAbiolinksGUI)
# TCGAbiolinksGUI()
de <- read.table("DE_FDR.xls", header = T)
de <- de[,c(1,2,6,7)]
dm <- read.csv("DM_nearest_target_position_value.csv", header = T)
dm <- dm[,c(3,4,5,6,7,8)]
colnames(dm)[4] <- 'gene_id'
dm_island <- dm[dm$Feature_Type == "Island",]

de_dm_island <- merge(x = de, y = dm_island, 'gene_id', y.all = T)
library(dplyr)
de_dm_island <- de_dm_island %>%
  mutate(Gene = 
           case_when(log2FoldChange > 1 & mean.Primary.Tumor.minus.mean.Solid.Tissue.Normal < -0.25 ~ 'G+M-',
                     log2FoldChange < -1 & mean.Primary.Tumor.minus.mean.Solid.Tissue.Normal > 0.25 ~ 'G-M+',))
de_dm_island$Gene[is.na(de_dm_island$Gene)] <- 'Others'
write.csv(de_dm_island, file = 'methylation-transcripts.csv')
library(ggplot2)
# set.seed(1234)
# ss <- sample(1:32, 15)
df <- de_dm_island

p <- ggplot(df, aes(x=mean.Primary.Tumor.minus.mean.Solid.Tissue.Normal, y=log2FoldChange, color = Gene, shape = position)) + 
     geom_point() +
     geom_vline(xintercept = c(-0.25,0.25), linetype = "dashed", color = "blue") +
     geom_hline(yintercept = c(-1,1), linetype="dashed", color = "blue") + 
     labs(x = "CpG Methylation (Delta of Mean Beta Value)", y = "Gene Expression (Log2FoldChange)") +
     scale_color_manual(values = c("G-M+" = "orange", "G+M-" = "firebrick2", "Others" = "gray45"))


require("ggrepel")
# set.seed(42)
p + geom_text_repel(data = subset(df, Gene == "G+M-" | Gene == "G-M+"), aes(label=gene_id), size = 3.5, show.legend=F)





# plot(data = de_dm_island, mean.Primary.Tumor.minus.mean.Solid.Tissue.Normal, log2FoldChange, main="Scatterplot Example",
#      xlab="mean of DM CpG ", ylab="log2FC of DE gene ", pch=19)


