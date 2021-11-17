########################################################################################
setwd("/Users/namuhanz/cancer/transcripts/compare DE/file")
### input tables and keep gene_ID and up-down columns #####################
c1 <- read.table("./c1.xls")
c4 <- read.table("./c4.xls")
c1 <- c1[-1,]
c1 <- c1[,c(1,7)]
colnames(c1) <- c("gene_id", "sig")
c4 <- c4[-1,]
c4 <- c4[,c(1,7)]
colnames(c4) <- c("gene_id", "sig")
### check the overlap between DE list 1 and DE list 4 #####################
num_gene_in_c1 <- nrow(c1)
num_gene_in_c4 <- nrow(c4)

c1_c4 <- intersect(c1$gene_id,c4$gene_id)
num_overlap <- length(c1_c4)

c1_c4 <- data.frame(c1_c4)

write.csv(c1_c4, file = "./c1-c4 overlap.csv")
print(c(num_gene_in_c1, num_gene_in_c4, num_overlap))



