setwd("/Users/namuhanz/cancer/miRNA/DE comparison/file")
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
female_specific_up <- female_specific[female_specific$sig=="Up",]
female_specific_down <- female_specific[female_specific$sig=="Down",]
write.csv(female_specific_up[,c(1,3)], file = "./female_specific_up.csv", row.names = FALSE)
write.csv(female_specific_down[,c(1,3)], file = "./female_specific_down.csv", row.names = FALSE)
write.csv(female_specific[,c(1,3)], file = "./female_specific.csv", row.names = FALSE)

male_specific <- com[com$chng_type=="-",]
male_specific_up <- male_specific[male_specific$sig=="Up",]
male_specific_down <- male_specific[male_specific$sig=="Down",]
write.csv(male_specific_up[,c(1,3)], file = "./male_specific_up.csv", row.names = FALSE)
write.csv(male_specific_down[,c(1,3)], file = "./male_specific_down.csv", row.names = FALSE)
write.csv(male_specific[,c(1,3)], file = "./male_specific.csv", row.names = FALSE)

print(c(
  nrow(female_specific_up),
  nrow(female_specific_down),
  nrow(male_specific_up),
  nrow(male_specific_down)
))
