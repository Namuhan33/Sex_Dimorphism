setwd("/Users/namuhanz/cancer/transcripts/enrichment analysis/file")
library(clusterProfiler)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)

genelist <- list(female_bladder=data.frame(),male_bladder=data.frame(),
                 female_colon=data.frame(),male_colon=data.frame(),
                 female_kidney=data.frame(), male_kidney=data.frame(),
                 female_liver=data.frame(),male_liver=data.frame(),
                 female_stomach=data.frame(),male_stomach=data.frame(),
                 female_thyroid=data.frame(),male_thyroid=data.frame())

for (i in 1:length(genelist)) {
  genelist[[i]] <- read.csv(paste("./",names(genelist[i]),".csv",sep = ""), header = T) 
}

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

for(j in 1:length(genelist)){
  symbol <- genelist[[j]]$gene_id
  entrezid <- mapIds(x = org.Hs.eg.db,
                     keys = symbol,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
  em <- enricher(entrezid, TERM2GENE=m_t2g)
  msigdb <- em@result
  path <- paste(names(genelist[j]), "_hallmarks.csv",sep = "")
  write.csv(msigdb, file = path)
}

  
  
  
  
  
