setwd("/Users/namuhanz/cancer/methylation/CG classification")
library(biomaRt)
# library(org.Hs.eg.db)

fulllist <- read.csv("jhu-usc.edu_KICH.HumanMethylation450.1.lvl-3.TCGA-KO-8415-01A-11D-2312-05.gdc_hg38.txt", sep = "\t", header = T)
# fulllist <- fulllist[fulllist$Gene_Symbol != ".",]

library(FDb.InfiniumMethylation.hg19)
hm450 = get450k()
cpg = fulllist$Composite.Element.REF
probes = hm450[cpg]
nearest_gene <- getNearestGene(probes)
nearest_TSS <-  getNearestTSS(probes)
nearest_transcript <- getNearestTranscript(probes)

# cgconvert <- read.csv("./cgconvert.csv", header = T)
# cgconvert <- cgconvert[,-1]
# cgconvert <- cgconvert[cgconvert$gene_symbols != '.',]
# cgconvert <- cbind(cgconvert, fulllist[,c("Start","End")])
# 
# index_multipletar <-  grep("., .", cgconvert$gene_symbols, ignore.case = FALSE, perl = FALSE, value = FALSE,
#                             fixed = FALSE, useBytes = FALSE, invert = FALSE)
# rownames(cgconvert) <- NULL
# multipletar <- cgconvert[index_multipletar,]
# 
# cpg <- read.csv("./cpgislandposition.csv", header = T)
# cpg <- cpg[,-1]
# cpg <- cpg[-1,]
# genesymbols <- unique(cpg$Gene_Symbol)





ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# filters <- listFilters(ensembl)
filters_1 <- 'external_gene_name'
# attributes <- listAttributes(ensembl)
# attributes_1 <- 'entrezgene'
# attributes_2 <- 'start_position'
# attributes_3 <- 'chromosome_name'
# 
symbol_to_transcript_length <- getBM(attributes = c(filters_1, "ensembl_transcript_id", "transcript_length", "chromosome_name","transcription_start_site", "5_utr_end"),
                                     filters = filters_1,
                                     values = unique(nearest_gene$nearestGeneSymbol),
                                     mart = ensembl)

write.csv(symbol_to_transcript_length, file = 'all transcripts.csv')
# 
# symbol_to_startpoint <- getBM(attributes = c(filters_1, attributes_2, 'strand', 'chromosome_name'),
#                               filters = filters_1,
#                               values = genesymbols, 
#                               mart = ensembl)
symbol_to_transcript_length <- read.csv('all transcripts.csv', header = T)
symbol_to_transcript_length <- symbol_to_transcript_length[,-1]
library(dplyr)
# s_t_t_l <- filter(symbol_to_transcript_length, "chromosome_name" != "CHR_^")
# 
# vars <- "chromosome_name"
# s_t_t_l<- symbol_to_transcript_length %>%
#   filter(
#     .data[[vars]] != 'CHR_^',
#   )

s_t_t_l <- filter(symbol_to_transcript_length, substr(chromosome_name,1,3) != "CHR")
nearest_gene <- cbind(rownames(nearest_gene), data.frame(nearest_gene, row.names=NULL))

require(data.table)
group <- as.data.table(s_t_t_l)

result <- group[group[, .I[which.max(transcript_length)], by=external_gene_name]$V1]

result_1 <- merge(nearest_gene, result, by.x="nearestGeneSymbol", by.y="external_gene_name", all.x=TRUE)

result_1 <- result_1[,c(-3,-4,-5)]
colnames(result_1)[c(1,2)] <- c("genesymbol", "cgID")
result_2 <- result_1[order(result_1$cgID), ]
rownames(result_2) <- NULL

result_3 <- cbind(result_2, fulllist[,c(3,4,5,10,11)])
result_3 <- result_3[,c(2,1,3:12)]

result_3$Chromosome<-gsub("chr","",as.character(result_3$Chromosome))
result_3$posi <- (result_3$Start - result_3$transcription_start_site)
result_4 <- result_3[result_3$posi >= -1500,]
result_5 <- result_4[result_4$chromosome_name == result_4$Chromosome,]
library(tidyr)
result_5 <- result_5 %>% drop_na(cgID)

write.csv(result_5, file = "reference.csv")
# InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
# met <- metadata(FDb.InfiniumMethylation.hg19) ## need to fetch genome
# genome(InfiniumMethylation) <- met[which(met[,'name']=='Genome'),'value']
# InfiniumMethylation <- sort(InfiniumMethylation)
# 
# data(hg19.islands)
# CGI.probes <- subsetByOverlaps(InfiniumMethylation, hg19.islands)
# head(CGI.probes)
# tail(CGI.probes)

ref <- read.csv("reference.csv", header = T)
ref <- ref[,-1]

ref <- ref %>% 
  mutate(position = if_else(posi >= 0, 'E', 'P'))

ref <- ref[,c(1,2,12,14)]
write.csv(ref, file = "ref.csv")




















