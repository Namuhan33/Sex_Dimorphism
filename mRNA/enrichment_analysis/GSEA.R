### tutorial: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
setwd("/Users/namuhanz/cancer/transcripts/enrichment analysis/file")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)
library(edgeR)

d_f <- read.csv("female_genelist.csv")
d_m <- read.csv("male_genelist.csv")

### feature 1: numeric vector
geneList_f = d_f[,3]
geneList_m = d_m[,3]
### feature 2: named vector
names(geneList_f) = as.character(d_f[,2])
names(geneList_m) = as.character(d_m[,2])
### feature 3: decreasing orde
geneList_f = sort(geneList_f, decreasing = TRUE)
geneList_m = sort(geneList_m, decreasing = TRUE)

gse_f <- gseGO(geneList=geneList_f, 
             ont ="ALL", 
             # nPerm = 10000,
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             pAdjustMethod = "BH",
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db,
             eps = 0)
             #eps = 0

gse_m <- gseGO(geneList=geneList_m, 
               ont ="ALL", 
               # nPerm = 10000,
               keyType = "SYMBOL", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               pAdjustMethod = "BH",
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db,
               eps = 0)
               # eps = 0

dotplot(gse_f, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot(gse_m, showCategory=10, split=".sign") + facet_grid(.~.sign)

xf <- pairwise_termsim(gse_f)
xm <- pairwise_termsim(gse_m)
par(mar = c(5.1, 4.1, 4.1, 2.1))
emapplot(xf, showCategory = 10)
emapplot(xm, showCategory = 10)

cnetplot(gse_f, categorySize="pvalue", foldChange=geneList_f, showCategory = 3)
cnetplot(gse_m, categorySize="pvalue", foldChange=geneList_m, showCategory = 3)

ridgeplot(gse_f) + labs(x = "enrichment distribution")
par(mar = c(100, 4.1, 4.1, 100))
ridgeplot(gse_m) + labs(x = "enrichment distribution")

write.csv(gse_f, file = "./GSEA_female_specific.csv")
write.csv(gse_m, file = "./GSEA_male_specific.csv")




