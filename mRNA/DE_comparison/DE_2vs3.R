########################################################################################
setwd("/Users/namuhanz/cancer/transcripts/compare DE/file")
library(compareDF)
### input tables and keep gene_ID and up-down columns #####################
c2 <- read.table("./c2.xls")
c3 <- read.table("./c3.xls")
c2 <- c2[-1,]
c2 <- c2[,c(1,7)]
colnames(c2) <- c("gene_id", "sig")
c3 <- c3[-1,]
c3 <- c3[,c(1,7)]
colnames(c3) <- c("gene_id", "sig")

com = compare_df(c2,c3,c("gene_id"))$comparison_df
female_specific <- com[com$chng_type=="+",]
write.csv(female_specific[,c(1,3)], file = "./female_specific.csv", row.names = FALSE)
male_specific <- com[com$chng_type=="-",]
write.csv(male_specific[,c(1,3)], file = "./male_specific.csv", row.names = FALSE)

print(c(nrow(c2), nrow(c3), nrow(female_specific), nrow(male_specific)))
