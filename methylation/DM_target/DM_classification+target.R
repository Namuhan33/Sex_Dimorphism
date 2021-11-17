setwd("/Users/namuhanz/cancer/methylation/DM target/file")
dm <- read.csv("DMR_results_sample.type_Primary.Tumor_Solid.Tissue.Normal_pcut_0.05_meancut_0.25.csv", header = T)

cg_convert <- read.csv("../ref.csv", header = T)
cg_convert <- cg_convert[,-1]
dm <- dm[dm$status != 'Not Significant',]
dm <- dm[,c(1,4,6,7)]
dm$status[dm$status == 'Hypomethylated in Primary Tumor'] <- 'hypo'
dm$status[dm$status == 'Hypermethylated in Primary Tumor'] <- 'hyper'
colnames(dm)[1] <- 'cgID'
rownames(dm) <- NULL

print(nrow(dm[dm$status == 'hyper',]))
print(nrow(dm[dm$status == 'hypo',]))

library(tidyr)
df_dm <- merge(x = dm, y = cg_convert,'cgID', all.x = T)
df_dm <- df_dm %>% drop_na()

df_dm$Feature_Type[df_dm$Feature_Type != "Island"] <- 'OffIsland'

write.csv(df_dm, file = './DM_nearest_target_position_value.csv')

bp <- df_dm[,c(1,4,6,7)]
rownames(bp) <- bp[,1]
bp <- bp[,-1]
res1 <- bp[,c(2,1,3)]

I_H <- res1[res1$Feature_Type == "Island" & res1$status == "hyper",]
I_H_E <- nrow(I_H[I_H$position == "E",]) 
I_H_P <- nrow(I_H[I_H$position == "P",]) 

O_H <- res1[res1$Feature_Type == "OffIsland" & res1$status == "hyper",]
O_H_E <- nrow(O_H[O_H$position == "E",]) 
O_H_P <- nrow(O_H[O_H$position == "P",])

I_h <- res1[res1$Feature_Type == "Island" & res1$status == "hypo",]
I_h_E <- nrow(I_h[I_h$position == "E",])
I_h_P <- nrow(I_h[I_h$position == "P",])

O_h <- res1[res1$Feature_Type == "OffIsland" & res1$status == "hypo",]
O_h_E <- nrow(O_h[O_h$position == "E",])
O_h_P <- nrow(O_h[O_h$position == "P",])

TEST <- data.frame('Status' = c('hyper','hyper','hyper','hyper','hypo','hypo','hypo','hypo'), 
                    'Location'=c('E','E','P','P','E','E','P','P'), 
                    'Feature' = rep(c('Island', 'Outside'),4))
TEST$value <- c(I_H_E, O_H_E, I_H_P, O_H_P, I_h_E, O_h_E, I_h_P, O_h_P)

library(ggplot2)
ggplot(data=TEST, aes(y = value, x = Status, fill = Location)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() +
  facet_grid( ~ Feature) +
  labs(title = "Classification of all mapped DM sites") +
  geom_text(aes(label = value), position=position_dodge(width=0.9), vjust=-0.25)
 
  

# test <- TEST[TEST$Feature == 'Island',]
# 
# ggplot() +
#   geom_bar(data=test, aes(y = value, x = Status, fill = Location), stat="identity",
#            position='stack') +
#   theme_bw() + 
#   labs(title = "Classification of all mapped DM sites in the CpG island")




# tt <- rbind(TEST,TEST)
# tt$type <- c(rep('male',8),rep('female',8))
# 
# ggplot(data=tt, aes(y = value, x = Status, fill = Location)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   theme_bw() + 
#   facet_grid( type~ Feature) +
#   labs(title = "Classification of all mapped DM sites") +
#   geom_text(aes(label = value), position=position_dodge(width=0.9), vjust=-0.25)














# 
# rownames(dm_hyper) <- NULL
# rownames(dm_hypo) <- NULL
# colnames(dm_hyper)[1] <- 'cgID'
# colnames(dm_hypo)[1] <- 'cgID'
# 
# library(tidyr)
# df_hyper <- merge(x = dm_hyper, y = cg_convert,'cgID', all.x = T)
# df_hyper <- df_hyper %>% drop_na()
# 
# df_hypo <- merge(x = dm_hypo, y = cg_convert,'cgID', all.x = T)
# df_hypo <- df_hypo %>% drop_na()
# 
# inside_island_hyper <- df_hyper[df_hyper$Feature_Type == "Island",]
# outside_island_hyper <- df_hyper[df_hyper$Feature_Type != "Island",]
# inside_island_hypo <- df_hypo[df_hypo$Feature_Type == "Island",]
# outside_island_hypo <- df_hypo[df_hypo$Feature_Type != "Island",]
# 
# exon_hyper_inisland <- inside_island_hyper[df_hyper$position == 'E', ]
# promoter_hyper_inisland <- inside_island_hyper[df_hyper$position == 'P', ]
# 
# exon_hypo_inisland <- outside_island_hyper[df_hypo$position == 'E', ]
# promoter_hypo_inisland <- outside_island_hyper[df_hypo$position == 'P', ]
# 
# exon_hyper_outisland <- outside_island_hyper[df_hyper$position == 'E', ]
# promoter_hyper_outisland <- outside_island_hyper[df_hyper$position == 'P', ]
# 
# exon_hypo_outisland <- outside_island_hyper[df_hypo$position == 'E', ]
# promoter_hypo_outisland <- outside_island_hyper[df_hypo$position == 'P', ]










