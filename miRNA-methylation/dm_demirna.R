setwd("/Users/namuhanz/cancer/miRNA-methylation/file")
ref <- read.csv("../ref_cg_mirna.csv", header = T)
ref <- ref[,-1]
ref <- ref[,c(1,5,6,7,8)]
colnames(ref) <- c('cgID', 'miRNA', 'accession','position', 'island')
mir <- read.table("./DEmiRNA_FDR.xls", header = T, sep = "")
dm <- read.csv("./DMR_results_sample.type_Primary.Tumor_Solid.Tissue.Normal_pcut_0.05_meancut_0.25.csv", header = T)
dm <- dm[!dm$status=='Not Significant',]
dm <- dm[,c(1,4,6,7)]
colnames(dm) <- c('cgID', 'delta.mean', 'p.adj', 'sig')
dm$sig[dm$sig=='Hypomethylated in Primary Tumor'] <- 'hypo'
dm$sig[dm$sig=='Hypermethylated in Primary Tumor'] <- 'hyper'

m1 <- merge(x = ref, y = dm, by = 'cgID', all.x = T)

library(tidyr)
m1 <- m1 %>% drop_na()

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      version = 'GRCh37')
m2 <- m1$accession
m3 <- getBM(attributes = c('refseq_ncrna', 'mirbase_id'),
            filters = 'refseq_ncrna',
            values = m2, 
            mart = ensembl)
m4 <- merge(x=m1, y=m3, by.x = 'accession', by.y = 'refseq_ncrna', all = T)
m5 <- merge(x=m4, y=mir, by.x = 'mirbase_id', by.y = 'gene_id', all = T)


m5 <- m5 %>% drop_na()

write.csv(m5, file = 'DMC_DEmiRNA.csv')



library(multiMiR)
library(miRNAmeConverter)
nc <- MiRNANameConverter()
new_id <- translateMiRNAName(nc, m5$mirbase_id, versions = c(17,20,21),sequenceFormat = 1,verbose = F)
mir_id <- unique(c(new_id$input,new_id$v17.0,new_id$v20.0,new_id$v21.0, m5$mirbase_id))
res <- get_multimir(org     = 'hsa',
                    mirna   = mir_id,
                    table   = 'validated',
                    summary = TRUE)
mir_target <- res@summary[, c(2,3)]
mir_target <- mir_target[mir_target$target_symbol != '',]
mir_target <- unique(mir_target)
miR_df <- new_id[,c(2,3,4,5)]
miR_df$input <- tolower(miR_df$input)
missing <- setdiff(m5$mirbase_id, miR_df$input)
df_temp <- data.frame(input = missing, v17.0 = missing, v20.0 = missing, v21.0 = missing)
miR_df <- rbind(miR_df, df_temp)
rownames(miR_df) <-NULL
df_1 <- miR_df[,c(1,2)]
df_2 <- miR_df[,c(1,3)]
df_3 <- miR_df[,c(1,4)]
m <- merge(x = mir_target, y = df_1, by.x = 'mature_mirna_id', by.y = 'v17.0', all.x = T)
m <- merge(x = m, y = df_2, by.x = 'mature_mirna_id', by.y = 'v20.0', all.x = T)
m <- merge(x = m, y = df_3, by.x = 'mature_mirna_id', by.y = 'v21.0', all.x = T)
library(dplyr)

m <- m %>% 
  mutate(miRNA = coalesce(input.x, input.y, input))

m <- m[,c(1,2,6)]

m6 <- merge(x=m, y=m5, by.x = 'miRNA', by.y = 'mirbase_id', all.y = T)
m6 <- m6[,c(1,3,7,8,9,10,11,12,17)]

deg <- read.csv("./DE_FDR.xls", sep = "", header = T)
m7  <- merge(x= deg, y=m6, by.x = 'gene_id', by.y = 'target_symbol', all.x = T)
m7 <- m7[,c(-3,-4,-5)]
m7 <- m7 %>% drop_na()

write.csv(m7, file = 'dm_demirna_de.csv')

colnames(m7) <- c('gene','l2fc.gene','p.adj.gene','sig.gene','mirna','position','island','delta.mean','p.adj.cpg','sig.cpg','l2fc.mirna','sig.mirna')

library("scatterplot3d")
m8 <- m7[,c(1,2,6,8,11)]
write.csv(m8,file = '3D scatter.csv')

shapes <- c(16,17,18)
shapes <- shapes[as.factor(m8$position)]
colors <- c("#999999", "#E69F00", "#56B4E9")
colors <- colors[as.factor(m8$position)]
s3d <- scatterplot3d(m8[,c(2,4,5)],
              angle = 55,
              main = 'DM-DEmiR-DEG',
              xlab = 'Gene Expression (log2FC)',
              ylab = 'Methylation (delta mean)',
              zlab = 'miRNA Expression (log2FC)',
              pch = 16, 
              color = colors
              )
legend('bottom', legend = levels(as.factor(m8$position)),
       col =  c("#999999", "#E69F00", "#56B4E9"), pch = 16,
       inset = -0.2, xpd = TRUE, horiz = TRUE)


text(s3d$xyz.convert(m8[,c(2,4,5)]),
     labels = m8$gene,
     cex = 0.5
     )









