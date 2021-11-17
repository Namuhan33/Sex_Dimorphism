setwd("/Users/namuhanz/cancer/mRNA-miRNA/file")
library(multiMiR)
library(miRNAmeConverter)
demir <- read.table("./DEmiRNA_FDR&FC.xls", header = T)
detranscript <- read.table("./DE_FDR.xls", header = T)
nc <- MiRNANameConverter()
new_id <- translateMiRNAName(nc, demir$gene_id, versions = c(17,20,21),sequenceFormat = 1,verbose = F)
mir_id <- unique(c(new_id$input,new_id$v17.0,new_id$v20.0,new_id$v21.0, demir$gene_id))
res <- get_multimir(org     = 'hsa',
                    mirna   = mir_id,
                    table   = 'validated',
                    summary = TRUE)

mir_target <- res@summary[, c(2,3)]
mir_target <- mir_target[mir_target$target_symbol != '',]
mir_target <- unique(mir_target)
miR_df <- new_id[,c(2,3,4,5)]

miR_df$input <- tolower(miR_df$input)
missing <- setdiff(demir$gene_id, miR_df$input)
df_temp <- data.frame(input = missing, v17.0 = missing, v20.0 = missing, v21.0 = missing)
miR_df <- rbind(miR_df, df_temp)
rownames(miR_df) <-NULL
df_1 <- miR_df[,c(1,2)]
df_2 <- miR_df[,c(1,3)]
df_3 <- miR_df[,c(1,4)]

m1 <- merge(x = mir_target, y = df_1, by.x = 'mature_mirna_id', by.y = 'v17.0', all.x = T)
m1 <- merge(x = m1, y = df_2, by.x = 'mature_mirna_id', by.y = 'v20.0', all.x = T)
m1 <- merge(x = m1, y = df_3, by.x = 'mature_mirna_id', by.y = 'v21.0', all.x = T)

library(dplyr)

m1 <- m1 %>% 
  mutate(miRNA = coalesce(input.x, input.y, input))

m1 <- m1[,c(1,2,6)]


m2 <- merge(x = m1, y = demir, by.x = 'miRNA', by.y = 'gene_id', all.x = T)

m3 <- merge(x = m2, y = detranscript, by.x = 'target_symbol', by.y = 'gene_id', all.y = T)

library(tidyr)
m3 <- m3 %>% drop_na()

m3 <- m3[,c(1,2,4,9,10,15)]
colnames(m3) <- c('gene','miRNA','l2fc_miRNA','sig_miRNA', 'l2fc_T', 'sig_T')


m3 <- m3 %>%
  mutate(correlation = 
           case_when(l2fc_T > 1 & l2fc_miRNA < -1 ~ 'G+miR-',
                     l2fc_T < -1 & l2fc_miRNA > 1 ~ 'G-miR+',))
m3$correlation[is.na(m3$correlation)] <- 'Others'
write.csv(m3, file = 'miRNA-transcripts.csv')

library(ggplot2)
# set.seed(1234)
# ss <- sample(1:32, 15)
df <- m3

p <- ggplot(df, aes(x=l2fc_miRNA, y=l2fc_T, color = correlation)) + 
  geom_point(shape = '.') +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = c(-1,1), linetype="dashed", color = "blue") + 
  labs(x = "miRNA Expression (Log2FoldChange)", y = "Gene Expression (Log2FoldChange)") +
  scale_color_manual(values = c("G-miR+" = "orange", "G+miR-" = "firebrick2", "Others" = "gray45"))


require("ggrepel")
p + geom_text_repel(data = subset(df, correlation == 'G+miR-' | correlation == 'G-miR+'), aes(label=gene), size = 2.5, show.legend=F)







