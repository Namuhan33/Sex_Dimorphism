setwd("/Users/namuhanz/cancer/mRNA-miRNA/file")
# Sex-specific DE genes #########
genelist <- list(female_bladder=data.frame(),male_bladder=data.frame(),
                 female_colon=data.frame(),male_colon=data.frame(),
                 female_kidney=data.frame(), male_kidney=data.frame(),
                 female_liver=data.frame(),male_liver=data.frame(),
                 female_stomach=data.frame(),male_stomach=data.frame(),
                 female_thyroid=data.frame(),male_thyroid=data.frame())

for (i in 1:length(genelist)) {
  genelist[[i]] <- read.csv(paste("./",names(genelist[i]),".csv",sep = ""), header = T) 
}

genes <- list(f_bladder=character(),m_bladder=character(),
              f_colon=character(),m_colon=character(),
              f_kidney=character(),m_kidney=character(),
              f_liver=character(),m_liver=character(),
              f_stomach=character(),m_stomach=character(),
              f_thyroid=character(),m_thyroid=character())
for (k in 1:length(genelist)) {
  genes[[k]] <- genelist[[k]]$gene_id
}

genes_f <- genes[c(1,3,5,7,9,11)]
genes_m <- genes[c(2,4,6,8,10,12)]

xxx_1 <- Reduce(intersect, list(genes_m[[1]], genes_m[[3]], genes_m[[4]], genes_m[[5]], genes_m[[6]]))
xxx_2 <- Reduce(intersect, list(genes_m[[1]], genes_m[[3]], genes_m[[2]], genes_m[[5]], genes_m[[6]]))
xxx_3 <- Reduce(intersect, list(genes_m[[1]], genes_m[[4]], genes_m[[2]], genes_m[[5]], genes_m[[6]]))
xxx_4 <- Reduce(intersect, list(genes_f[[4]], genes_f[[2]], genes_f[[3]], genes_f[[6]]))
xxx_5 <- Reduce(intersect, list(genes_f[[1]], genes_f[[2]], genes_f[[5]], genes_f[[6]]))
xxx_6 <- Reduce(intersect, list(genes_f[[3]], genes_f[[4]], genes_f[[5]], genes_f[[6]]))
# Sex-specific DE gene enrichment analysis against hallmark gene sets ###########
hallmarklist <- list(female_bladder_hallmarks=data.frame(),male_bladder_hallmarks=data.frame(),
                     female_colon_hallmarks=data.frame(),male_colon_hallmarks=data.frame(),
                     female_kidney_hallmarks=data.frame(), male_kidney_hallmarks=data.frame(),
                     female_liver_hallmarks=data.frame(),male_liver_hallmarks=data.frame(),
                     female_stomach_hallmarks=data.frame(),male_stomach_hallmarks=data.frame(),
                     female_thyroid_hallmarks=data.frame(),male_thyroid_hallmarks=data.frame())
for (j in 1:length(hallmarklist)) {
  hm <- read.csv(paste("./",names(hallmarklist[j]),".csv",sep = ""), header = T)
  hallmarklist[[j]] <- hm[hm$p.adjust <= 0.05,]
}

hallmarks <- list(f_bladder_h=character(), m_bladder_h=character(),
                  f_colon_h=character(), m_colon_h=character(),
                  f_kidney_h=character(), m_kidney_h=character(),
                  f_liver_h=character(), m_liver_h=character(),
                  f_stomach_h=character(), m_stomach_h=character(),
                  f_thyroid_h=character(), m_thyroid_h=character())

for (k in 1:length(hallmarklist)) {
  hallmarks[[k]] <- hallmarklist[[k]]$'ID'
}

hallmarks_f <- hallmarks[c(1,3,5,7,9,11)]
hallmarks_m <- hallmarks[c(2,4,6,8,10,12)]

xxx_7 <- Reduce(intersect, hallmarks_f[c(1,4)])
xxx_8 <- Reduce(intersect, hallmarks_m[c(2,3,4,6)])
xxx_9 <- Reduce(intersect, c(hallmarks_m[c(1,3,5)],hallmarks_f[c(1,4)]))

# Sex-specific DE miRNAs ############
miRNAlist <- list(female_bladder=data.frame(),male_bladder=data.frame(),
                  female_colon=data.frame(),male_colon=data.frame(),
                  female_kidney=data.frame(), male_kidney=data.frame(),
                  female_liver=data.frame(),male_liver=data.frame(),
                  female_stomach=data.frame(),male_stomach=data.frame(),
                  female_thyroid=data.frame(),male_thyroid=data.frame())
for (i in 1:length(miRNAlist)) {
  miRNAlist[[i]] <- read.csv(paste("./",names(miRNAlist[i]),".csv",sep = ""), header = T) 
}

miRNAs <- list(f_bladder_mi=character(), m_bladder_mi=character(),
               f_colon_mi=character(), m_colon_mi=character(),
               f_kidney_mi=character(), m_kidney_mi=character(),
               f_liver_mi=character(), m_liver_mi=character(),
               f_stomach_mi=character(), m_stomach_mi=character(),
               f_thyroid_mi=character(), m_thyroid_mi=character())

for (k in 1:length(miRNAlist)) {
  miRNAs[[k]] <- miRNAlist[[k]]$gene_id
}

miRNAs_f <- miRNAs[c(1,3,5,7,9,11)]
miRNAs_m <- miRNAs[c(2,4,6,8,10,12)]

xxx_10 <- Reduce(intersect, miRNAs_f[c(2,3,5)])
xxx_11 <- Reduce(intersect, miRNAs_f[c(1,2,4)])
xxx_12 <- Reduce(intersect, miRNAs_m[c(1,3,5)])
xxx_13 <- Reduce(intersect, miRNAs_m[c(2,3,6)])
xxx_14 <- Reduce(intersect, miRNAs_m[c(3,4,5)])
# Sex-specific DE miRNA true target ##########
truetargetlist <- list(female_bladder=data.frame(),male_bladder=data.frame(),
                       female_colon=data.frame(),male_colon=data.frame(),
                       female_kidney=data.frame(), male_kidney=data.frame(),
                       female_liver=data.frame(),male_liver=data.frame(),
                       female_stomach=data.frame(),male_stomach=data.frame(),
                       female_thyroid=data.frame(),male_thyroid=data.frame())
for (i in 1:length(truetargetlist)) {
  truetargetlist[[i]] <- read.csv(paste("./",names(truetargetlist[i])," true target",".csv",sep = ""), header = T) 
}

truetargets <- list(f_bladder_tt=character(), m_bladder_tt=character(),
                    f_colon_tt=character(), m_colon_tt=character(),
                    f_kidney_tt=character(), m_kidney_tt=character(),
                    f_liver_tt=character(), m_liver_tt=character(),
                    f_stomach_tt=character(), m_stomach_tt=character(),
                    f_thyroid_tt=character(), m_thyroid_tt=character())

for (k in 1:length(truetargetlist)) {
  truetargets[[k]] <- truetargetlist[[k]]$x
}

truetargets_f <- truetargets[c(1,3,5,7,9,11)]
truetargets_m <- truetargets[c(2,4,6,8,10,12)]

xxx_15 <- Reduce(intersect, truetargets_f[c(2,3,5)])
xxx_16 <- Reduce(intersect, truetargets_m[c(1,2,3,5,6)])
# Sex-specific DM ##########
dmlist <- list(female_bladder_DM=data.frame(),male_bladder_DM=data.frame(),
               female_colon_DM=data.frame(),male_colon_DM=data.frame(),
               female_kidney_DM=data.frame(), male_kidney_DM=data.frame(),
               female_liver_DM=data.frame(),male_liver_DM=data.frame(),
               female_thyroid_DM=data.frame(),male_thyroid_DM=data.frame())
for (i in 1:length(dmlist)) {
  dmlist[[i]] <- read.csv(paste("./",names(dmlist[i]),".csv",sep = ""), header = T)
}

dm <- list(f_bladder_dm=character(), m_bladder_dm=character(),
           f_colon_dm=character(), m_colon_dm=character(),
           f_kidney_dm=character(), m_kidney_dm=character(),
           f_liver_dm=character(), m_liver_dm=character(),
           f_thyroid_dm=character(), m_thyroid_dm=character())

for (k in 1:length(dmlist)) {
  dm[[k]] <- dmlist[[k]]$cgid
}

dm_f <- dm[c(1,3,5,7,9)]
dm_m <- dm[c(2,4,6,8,10)]

xxx_16 <- Reduce(intersect, dm_f[c(1,2,4)])
xxx_17 <- Reduce(intersect, dm_f[c(1,3,4)])
xxx_18 <- Reduce(intersect, dm_m[c(1,2,4,5)])
xxx_19 <- Reduce(intersect, dm_m[c(1,2,4,3)])
xxx_20 <- Reduce(intersect, dm_m[c(1,2,3,5)])
xxx_21 <- Reduce(intersect, dm_m[c(1,4,3,5)])
xxx_22 <- Reduce(intersect, dm_m[c(2,4,3,5)])

# sex-specific DM true targets ##########
dmttlist <- list(female_bladder_dmtt=data.frame(),male_bladder_dmtt=data.frame(),
                 female_colon_dmtt=data.frame(),male_colon_dmtt=data.frame(),
                 female_kidney_dmtt=data.frame(), male_kidney_dmtt=data.frame(),
                 female_liver_dmtt=data.frame(),male_liver_dmtt=data.frame(),
                 female_thyroid_dmtt=data.frame(),male_thyroid_dmtt=data.frame())
for (i in 1:length(dmttlist)) {
  dmttlist[[i]] <- read.csv(paste("./",names(dmttlist[i]),".csv",sep = ""), header = T)
}

dmtt <- list(f_bladder_dmtt=character(), m_bladder_dmtt=character(),
             f_colon_dmtt=character(), m_colon_dmtt=character(),
             f_kidney_dmtt=character(), m_kidney_dmtt=character(),
             f_liver_dmtt=character(), m_liver_dmtt=character(),
             f_thyroid_dmtt=character(), m_thyroid_dmtt=character())

for (k in 1:length(dmttlist)) {
  dmtt[[k]] <- dmttlist[[k]]$x
}

dmtt_f <- dmtt[c(1,3,5,7,9)]
dmtt_m <- dmtt[c(2,4,6,8,10)]

xxx_23 <- Reduce(intersect, dmtt_f[c(4,3)])
xxx_24 <- Reduce(intersect, dmtt_m[c(1,4,3)])



