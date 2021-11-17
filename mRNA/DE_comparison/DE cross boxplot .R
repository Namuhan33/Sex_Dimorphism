setwd("/Users/namuhanz/cancer/transcripts/compare DE/file")

library(edgeR)
library(ggplot2)

num_mn <- 10
num_mt <- 12
num_fn <- 9
num_ft <- 9

# filtering and normalization#############
malecount <- read.csv("./Counts_homo_male.csv", head= F)
femalecount <- read.csv("./Counts_homo_female.csv", head = F)
colnames(malecount) <- malecount[1,]
colnames(femalecount) <- femalecount[1,]
rownames(malecount) <- malecount[,1]
rownames(femalecount) <- femalecount[,1]
malecount <- malecount[-1,]
malecount <- malecount[,-1]
femalecount <- femalecount[-1,]
femalecount <- femalecount[,-1]
countsdata <- cbind(malecount,femalecount)
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)

sampledescription <- data.frame()
SampleID <- data.frame(colnames(countdata))
Sex <- data.frame(c(rep("male", num_mn), 
                    rep("male", num_mt), 
                    rep("female", num_fn), 
                    rep("female", num_ft)),stringsAsFactors = T)
Type <- data.frame(c(rep("normal", num_mn), 
                     rep("tumor", num_mt), 
                     rep("normal", num_fn), 
                     rep("tumor", num_ft)),stringsAsFactors = T)
sampledescription <- cbind(SampleID, Sex, Type)
colnames(sampledescription) <- c("SampleID", "Sex", "Type")
target <- sampledescription
rownames(target) <- sampledescription[,1]
target <- target[,-1]
target$Sex<- as.factor(target$Sex)
target$Type <- as.factor(target$Type)

group <- factor(paste0(sampledescription$Sex, ".", sampledescription$Type))
y <- DGEList(countdata, group = group)
design <- model.matrix(~Sex+Type+Sex:Type, data = target)
keep <- filterByExpr(y,design)
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

# exact the most significant gene in c5 and in sex-specific genes#########
c5 <- read.table("./DEG_FDR&FC_m.xls", header = T)
c5 <- c5[order(c5$FDR),]
female <- read.csv("./female-c5 overlap.csv", header = T)
male <- read.csv("./male-c5.csv", header = T)
c5_female <- data.frame(matrix(nrow = nrow(female),
                               ncol = ncol(c5), NA))
c5_male <- data.frame(matrix(nrow = nrow(male),
                             ncol = ncol(c5), NA))
for(i in 1:nrow(female)){
  gene <- female$x[i]
  c5_female[i,] <- c5[c5$gene_id==gene,]
}
colnames(c5_female) <- colnames(c5)

for(j in 1:nrow(male)){
  gene <- male$x[j]
  c5_male[j,] <- c5[c5$gene_id==gene,]
}
colnames(c5_male) <- colnames(c5)

c5_female <- c5_female[order(c5_female$FDR),]
c5_male <- c5_male[order(c5_male$FDR),]

c5_gene_1 <- c5$gene_id[1]
female_gene_1 <- c5_female$gene_id[1]
male_gene_1 <- c5_male$gene_id[1]

# significant gene expression level exact ###########
counts.cpm <- cpm(y, log = F)
counts <- log(counts.cpm + 1)
counts_c5 <- counts[c5_gene_1,]
counts_f <- counts[female_gene_1,]
counts_m <- counts[male_gene_1,]
sex <- c(rep("male", num_mn), 
         rep("male", num_mt), 
         rep("female", num_fn), 
         rep("female", num_ft))
type <- c(rep("normal", num_mn), 
          rep("tumor", num_mt), 
          rep("normal", num_fn), 
          rep("tumor", num_ft))
df_c5 <- data.frame(sex,type,counts_c5)
df_f <- data.frame(sex,type,counts_f)
df_m <- data.frame(sex,type,counts_m)
# boxplot ###########################################
female.normal.1 <- df_c5[df_c5$sex== 'female' & df_c5$type=='normal',]$counts_c5
female.tumor.1 <- df_c5[df_c5$sex== 'female' & df_c5$type=='tumor',]$counts_c5
male.normal.1 <- df_c5[df_c5$sex== 'male' & df_c5$type=='normal',]$counts_c5
male.tumor.1 <- df_c5[df_c5$sex== 'male' & df_c5$type=='tumor',]$counts_c5
list.1 <- list('female normal' = female.normal.1, 
               'female tumor' = female.tumor.1, 
               'male normal'= male.normal.1, 
               'male tumor' = male.tumor.1)
boxplot(list.1,
        main = paste('Expression of most significant gene\n', c5_gene_1, ' in C5', sep = ''),
        at = c(1,2,4,5),
        col = c("orange","red"),
        names = c("normal.F", "tumor.F", "normal.M", "tumor.M"),
        notch = TRUE)


female.normal.f <- df_f[df_f$sex== 'female' & df_f$type=='normal',]$counts_f
female.tumor.f <- df_f[df_f$sex== 'female' & df_f$type=='tumor',]$counts_f
male.normal.f <- df_f[df_f$sex== 'male' & df_f$type=='normal',]$counts_f
male.tumor.f <- df_f[df_f$sex== 'male' & df_f$type=='tumor',]$counts_f
list.f <- list('female normal' = female.normal.f, 
               'female tumor' = female.tumor.f, 
               'male normal'= male.normal.f, 
               'male tumor' = male.tumor.f)
boxplot(list.f,
        main = paste('Expression of most significant female-specific gene\n', female_gene_1, ' in C5', sep = ''),
        at = c(1,2,4,5),
        col = c("orange","red"),
        names = c("normal.F", "tumor.F", "normal.M", "tumor.M"),
        notch = TRUE)

female.normal.m <- df_m[df_m$sex== 'female' & df_m$type=='normal',]$counts_m
female.tumor.m <- df_m[df_m$sex== 'female' & df_m$type=='tumor',]$counts_m
male.normal.m <- df_m[df_m$sex== 'male' & df_m$type=='normal',]$counts_m
male.tumor.m <- df_m[df_m$sex== 'male' & df_m$type=='tumor',]$counts_m
list.m <- list('female normal' = female.normal.m, 
               'female tumor' = female.tumor.m, 
               'male normal'= male.normal.m, 
               'male tumor' = male.tumor.m)
boxplot(list.m,
        main = paste('Expression of most significant male-specific gene\n', male_gene_1, ' in C5', sep = ''),
        at = c(1,2,4,5),
        col = c("orange","red"),
        names = c("normal.F", "tumor.F", "normal.M", "tumor.M"),
        notch = TRUE)






# point plot ##################
# library(dplyr)
# aaa <- df_c5 %>% group_by(type,sex) %>% 
#   select(counts_c5) %>% 
#   summarise_all(list(mean = mean,
#                      std = sd, 
#                      min = min, 
#                      max = max))
# 
# 
# bbb <- df_f %>% group_by(type,sex) %>% 
#   select(counts_f) %>% 
#   summarise_all(list(mean = mean,
#                      std = sd, 
#                      min = min, 
#                      max = max))
# 
# ccc <- df_m %>% group_by(type,sex) %>% 
#   select(counts_m) %>% 
#   summarise_all(list(mean = mean,
#                      std = sd, 
#                      min = min, 
#                      max = max))
# 
# 
# 
# 
# 
# 
#  p_1 <-ggplot(aaa, aes(x=type, y=mean, fill=sex, group=sex)) + 
#       geom_point(aes(shape=sex, color=sex, size = 3))+
#       geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2,
#                 position=position_dodge(0.05))+
#       xlab("Sample type")+
#       ylab("Expression of the most significant gene in c5")+
#       geom_text(aes(label = c5_gene_1))
# 
# p_2 <-ggplot(bbb, aes(x=type, y=mean, fill=sex, group=sex)) + 
#   geom_point(aes(shape=sex, color=sex, size = 3))+
#   geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2,
#                 position=position_dodge(0.05))+
#   xlab("Sample type")+
#   ylab("Expression of the most significant gene (female)")+
#   geom_text(aes(label = female_gene_1))
# 
# p_3 <-ggplot(ccc, aes(x=type, y=mean, fill=sex, group=sex)) + 
#   geom_point(aes(shape=sex, color=sex), size = 3)+
#   geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2,
#                 position=position_dodge(0.05))+
#   xlab("Sample type")+
#   ylab("Expression of the most significant gene (male)")+
#   geom_text(aes(label = male_gene_1))







