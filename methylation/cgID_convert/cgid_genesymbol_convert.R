setwd("/Users/namuhanz/cancer/methylation/cgID convert")

cg <- read.delim("./jhu-usc.edu_THCA.HumanMethylation450.2.lvl-3.TCGA-EM-A1CV-11A-11D-A13Z-05.gdc_hg38.txt", 
                 header = TRUE, 
                 sep = "\t", 
                 quote = "\"",
                 dec = ".", 
                 fill = TRUE, comment.char = "")
cg_gene <- cg[,c(1,6)]
cg_gene <- cbind(cg_gene, NA)
colnames(cg_gene) <- c("cg_id", "symbols", "gene_symbols")

for(i in 1:nrow(cg_gene)){
  cg_gene[i,3] <- toString(unique(unlist(strsplit(cg_gene[i,2], ";"))))
}

cgconvert <- cg_gene[,c(1,3)]
write.csv(cgconvert, file = "./cgconvert.csv")



# rownames(cg_gene) <- cg[,1]
# cg_gene <- cg_gene[,-1]
# cg_gene <- cg_gene[,-1]






