xxx <- read.csv("/Users/namuhanz/Downloads/HumanMethylation450_15017482_v1-2.csv", header = F)
xxx <- xxx[c(-1,-2,-3,-4,-5,-6,-7),]
colnames(xxx) <- xxx[1,]
xxx <- xxx[-1,]
xxx <- xxx[,c(1,12,16,17,22,23,24,26)]
xxx <- xxx[!xxx$UCSC_RefGene_Name=="",]

library(dplyr)
library(stringr)
xxx <- xxx %>%
  filter(str_detect(UCSC_RefGene_Name, "MIR"))


library(splitstackshape)
yyy <- cSplit(
  xxx, c("UCSC_RefGene_Name", "UCSC_RefGene_Group", "UCSC_RefGene_Accession"),
  sep = ";", direction = "long", type.convert = F)

yyy <- yyy %>%
  filter(str_detect(UCSC_RefGene_Name, "MIR"))

write.csv(yyy, file = "ref_cg_mirna.csv")







