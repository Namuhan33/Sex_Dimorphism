setwd("/Users/namuhanz/cancer/transcripts/enrichment analysis/file")
library(clusterProfiler)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
female <- read.csv("./female_specific.csv", header = T)
male <- read.csv("./male_specific.csv",header = T)
# m_df <- msigdbr(species = "Homo sapiens")

# over-representation analysis ######
# whether the DE genes are involved in a biological process
female_symbol <- female$gene_id
female_entrezid <- mapIds(x = org.Hs.eg.db,
                        keys = female_symbol,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")

male_symbol <- male$gene_id
male_entrezid <- mapIds(x = org.Hs.eg.db,
                          keys = male_symbol,
                          column = "ENTREZID",
                          keytype = "SYMBOL",
                          multiVals = "first")

# m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
#   dplyr::select(gs_name, entrez_gene)
# em <- enricher(gene_entrezid, TERM2GENE=m_t2g)

m_t2g <- msigdbr(species = "Homo sapiens", category = 3) %>% 
  dplyr::select(gs_name, entrez_gene)

em_female <- enricher(female_entrezid, TERM2GENE=m_t2g)
em_male <- enricher(male_entrezid,TERM2GENE=m_t2g)

female_msigdb <- em_female@result
male_msigdb <- em_male@result

write.csv(female_msigdb, file = "./female_msigdb.csv")
write.csv(male_msigdb, file = "./male_msigdb.csv")










