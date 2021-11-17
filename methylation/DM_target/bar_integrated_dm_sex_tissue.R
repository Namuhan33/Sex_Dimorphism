# generate dataframe to store the classification of counts of dm #############
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


library(tidyr)
df_dm <- merge(x = dm, y = cg_convert,'cgID', all.x = T)
df_dm <- df_dm %>% drop_na()

df_dm$Feature_Type[df_dm$Feature_Type != "Island"] <- 'OffIsland'
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
TEST$type <- c(rep('bladder', 8))
TEST$sex <- c(rep('female',8))

write.csv(TEST, file = 'female.bladder.csv')






# integrated bar plot ########
library(tidyverse)
f.b <- read.csv('./female.bladder.csv', header = T)
m.b <- read.csv('./male.bladder.csv', header = T)
f.c <- read.csv('./female.colon.csv', header = T)
m.c <- read.csv('./male.colon.csv', header = T)
f.k <- read.csv('./female.kidney.csv', header = T)
m.k <- read.csv('./male.kidney.csv', header = T)
f.l <- read.csv('./female.liver.csv', header = T)
m.l <- read.csv('./male.liver.csv', header = T)
f.t <- read.csv('./female.thyroid.csv', header = T)
m.t <- read.csv('./male.thyroid.csv', header = T)

m2 <- bind_rows(f.b, m.b, f.c, m.c, f.k, m.k, f.l, m.l, f.t, m.t)
m2 <- m2[,-1]
write.csv(m2, file = 'integrated_dm.csv') 
# require("ggrepel")

ggplot(data=m2, aes(y = value, x = Status, fill = Location)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw() +
    facet_grid( type ~ sex + Feature) +
    labs(title = "Classification of all mapped DM sites") +
    geom_text(aes(label = value), position=position_dodge(width=0.9), vjust=0.2, size=3.5)





