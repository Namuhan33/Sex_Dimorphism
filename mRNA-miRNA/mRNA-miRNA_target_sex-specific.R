### sex-specific genes vs sex-specific miRNA targets
setwd("/Users/namuhanz/cancer/mRNA-miRNA/file")

female_gene <- read.csv("./female_specific.csv", header = T)
male_gene <- read.csv("./male_specific.csv", header = T)
female_gene_up <- female_gene[female_gene$sig=="Up",]$gene_id
female_gene_down <- female_gene[female_gene$sig=="Down",]$gene_id
male_gene_up <- male_gene[male_gene$sig=="Up",]$gene_id
male_gene_down <- male_gene[male_gene$sig=="Down",]$gene_id

female_miR_up <- read.csv("./target_femaleup.csv", header = T)
female_miR_down <- read.csv("./target_femaledown.csv", header = T)
male_miR_up <- read.csv("./target_maleup.csv", header = T)
male_miR_down <- read.csv("./target_maledown.csv", header = T)

female_target_up <- female_miR_up$x
female_target_down <- female_miR_down$x
male_target_up <- male_miR_up$x
male_target_down <- male_miR_down$x

### say as genes
female_up <- intersect(female_gene_up, female_target_down)
female_down <- intersect(female_gene_down, female_target_up)
male_up <- intersect(male_gene_up,male_target_down)
male_down <- intersect(male_gene_down, male_target_up)

write.csv(female_up, file = "./female up.csv")
write.csv(female_down, file = "./female down.csv")
write.csv(male_up, file = "./male up.csv")
write.csv(male_down, file = "./male down.csv")

print(c(length(female_up), length(female_down),length(male_up),length(male_down)))




