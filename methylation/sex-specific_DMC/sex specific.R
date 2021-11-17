setwd("/Users/namuhanz/cancer/methylation/DM comparison 2vs3/file")
library(compareDF)

c2 <- read.csv("./c2.csv", header = T)
c3 <- read.csv("./c3.csv", header = T)

c2 <- c2[,-1]
c3 <- c3[,-1]

com = compare_df(c2,c3,c("cgid"))$comparison_df
female_specific <- com[com$chng_type=="+",]
write.csv(female_specific[,c(1,3)], file = "./female_specific_DM.csv", row.names = FALSE)
male_specific <- com[com$chng_type=="-",]
write.csv(male_specific[,c(1,3)], file = "./male_specific_DM.csv", row.names = FALSE)

print(c(nrow(c2), nrow(c3), nrow(female_specific), nrow(male_specific)))






