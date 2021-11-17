########################################################################################
setwd("/Users/namuhanz/cancer/transcripts/compare DE/file")
### input tables and keep gene_ID and up-down columns #####################
c1 <- read.table("./c1.xls")
c6 <- read.table("./c6.xls")
c1 <- c1[-1,]
c1 <- c1[,c(1,7)]
colnames(c1) <- c("gene_id", "sig")
c6 <- c6[-1,]
c6 <- c6[,c(1,7)]
colnames(c6) <- c("gene_id", "sig")
### check the overlap between DE list 1 and DE list 4 #####################
num_gene_in_c1 <- nrow(c1)
num_gene_in_c6 <- nrow(c6)

c1_c6 <- intersect(c1$gene_id,c6$gene_id)
num_overlap <- length(c1_c6)

c1_c6 <- data.frame(c1_c6)

write.csv(c1_c6, file = "./c1-c6 overlap.csv")
print(c(num_gene_in_c1, num_gene_in_c6, num_overlap))



