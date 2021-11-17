setwd("/Users/namuhanz/cancer/transcripts/compare DE/file")

c1_df <- read.table("./c1.xls", header = T)
c5_df <- read.table("./c5.xls", header = T)
female_df <- read.csv("./female_specific.csv", header = T)
male_df <- read.csv("./male_specific.csv", header = T)

c1 <- c1_df$gene_id
c5 <- c5_df$gene_id
female <- female_df$gene_id
male <- male_df$gene_id

c1_c5 <- intersect(c1,c5)
female_c5 <- intersect(female,c5)
male_c5 <- intersect(male,c5)

write.csv(c1_c5, file = "c1-c5 overlap.csv")
write.csv(female_c5, file = "female-c5 overlap.csv")
write.csv(male_c5, file = "./male-c5.csv")

print(c(length(c1), length(c1_c5), 
        length(female), length(female_c5), 
        length(male), length(male_c5)))







