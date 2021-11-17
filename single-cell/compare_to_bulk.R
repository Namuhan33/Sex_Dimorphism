setwd("/Users/namuhanz/cancer/single cell")
female_epi <- read.csv('female_Epithelial.csv', header = T)
male_epi <- read.csv('male_Epithelial.csv', header = T)
female_sc <- read.csv('female_all_de.csv', header = T)
male_sc <- read.csv('male_all_de.csv', header = T)
female_bulk <- read.csv('c2.xls', header = T, sep = '\t')
male_bulk <- read.csv('c3.xls', header = T, sep = '\t')
female_epi <- female_epi[abs(female_epi$avg_log2FC)>=1,]
female_sc <- female_sc[abs(female_sc$avg_log2FC)>=1,]
male_epi <- male_epi[abs(male_epi$avg_log2FC)>=1,]
male_sc <- male_sc[abs(male_sc$avg_log2FC)>=1,]
female_epi_bulk <- intersect(female_bulk$gene_id, female_epi$X)
female_sc_bulk <- intersect(female_bulk$gene_id, female_sc$X)
male_epi_bulk <- intersect(male_bulk$gene_id, male_epi$X)
male_sc_bulk <- intersect(male_bulk$gene_id, male_sc$X)
