setwd("/Users/namuhanz/cancer/mRNA-miRNA/file")
library(UpSetR)
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
# Sex-specific DM ##########
dmlist <- list(female_bladder_DM=data.frame(),male_bladder_DM=data.frame(),
               female_colon_DM=data.frame(),male_colon_DM=data.frame(),
               female_kidney_DM=data.frame(), male_kidney_DM=data.frame(),
               female_liver_DM=data.frame(),male_liver_DM=data.frame(),
               female_thyroid_DM=data.frame(),male_thyroid_DM=data.frame())
for (i in 1:length(dmlist)) {
  dmlist[[i]] <- read.csv(paste("./",names(dmlist[i]),".csv",sep = ""), header = T)
}
# sex-specific DM true targets ##########
dmttlist <- list(female_bladder_dmtt=data.frame(),male_bladder_dmtt=data.frame(),
               female_colon_dmtt=data.frame(),male_colon_dmtt=data.frame(),
               female_kidney_dmtt=data.frame(), male_kidney_dmtt=data.frame(),
               female_liver_dmtt=data.frame(),male_liver_dmtt=data.frame(),
               female_thyroid_dmtt=data.frame(),male_thyroid_dmtt=data.frame())
for (i in 1:length(dmttlist)) {
  dmttlist[[i]] <- read.csv(paste("./",names(dmttlist[i]),".csv",sep = ""), header = T)
}
# Upset plot ############
genes <- list(f_bladder=character(),m_bladder=character(),
              f_colon=character(),m_colon=character(),
              f_kidney=character(),m_kidney=character(),
              f_liver=character(),m_liver=character(),
              f_stomach=character(),m_stomach=character(),
              f_thyroid=character(),m_thyroid=character())

hallmarks <- list(f_bladder_h=character(), m_bladder_h=character(),
                  f_colon_h=character(), m_colon_h=character(),
                  f_kidney_h=character(), m_kidney_h=character(),
                  f_liver_h=character(), m_liver_h=character(),
                  f_stomach_h=character(), m_stomach_h=character(),
                  f_thyroid_h=character(), m_thyroid_h=character())

miRNAs <- list(f_bladder_mi=character(), m_bladder_mi=character(),
               f_colon_mi=character(), m_colon_mi=character(),
               f_kidney_mi=character(), m_kidney_mi=character(),
               f_liver_mi=character(), m_liver_mi=character(),
               f_stomach_mi=character(), m_stomach_mi=character(),
               f_thyroid_mi=character(), m_thyroid_mi=character())

truetargets <- list(f_bladder_tt=character(), m_bladder_tt=character(),
                    f_colon_tt=character(), m_colon_tt=character(),
                    f_kidney_tt=character(), m_kidney_tt=character(),
                    f_liver_tt=character(), m_liver_tt=character(),
                    f_stomach_tt=character(), m_stomach_tt=character(),
                    f_thyroid_tt=character(), m_thyroid_tt=character())

dm <- list(f_bladder_dm=character(), m_bladder_dm=character(),
           f_colon_dm=character(), m_colon_dm=character(),
           f_kidney_dm=character(), m_kidney_dm=character(),
           f_liver_dm=character(), m_liver_dm=character(),
           f_thyroid_dm=character(), m_thyroid_dm=character())

dmtt <- list(f_bladder_dmtt=character(), m_bladder_dmtt=character(),
             f_colon_dmtt=character(), m_colon_dmtt=character(),
             f_kidney_dmtt=character(), m_kidney_dmtt=character(),
             f_liver_dmtt=character(), m_liver_dmtt=character(),
             f_thyroid_dmtt=character(), m_thyroid_dmtt=character())

for (k in 1:length(dmlist)) {
  # genes[[k]] <- genelist[[k]]$gene_id
  # hallmarks[[k]] <- hallmarklist[[k]]$'ID'
  # miRNAs[[k]] <- miRNAlist[[k]]$gene_id
  # truetargets[[k]] <- truetargetlist[[k]]$x
  dm[[k]] <- dmlist[[k]]$cgid
  # dmtt[[k]] <- dmttlist[[k]]$x
}

genes_f <- genes[c(1,3,5,7,9,11)]
genes_m <- genes[c(2,4,6,8,10,12)]
hallmarks_f <- hallmarks[c(1,3,5,7,9,11)]
hallmarks_m <- hallmarks[c(2,4,6,8,10,12)]
miRNAs_f <- miRNAs[c(1,3,5,7,9,11)]
miRNAs_m <- miRNAs[c(2,4,6,8,10,12)]
truetargets_f <- truetargets[c(1,3,5,7,9,11)]
truetargets_m <- truetargets[c(2,4,6,8,10,12)]
dm_f <- dm[c(1,3,5,7,9)]
dm_m <- dm[c(2,4,6,8,10)]
dmtt_f <- dmtt[c(1,3,5,7,9)]
dmtt_m <- dmtt[c(2,4,6,8,10)]

library(ComplexHeatmap)
m1 <- make_comb_mat(genes, mode = "intersect")
m1 <- m1[comb_degree(m1) >= 2]
UpSet(m1)
ht = draw(UpSet(m1))
od = column_order(ht)
cs = comb_size(m1)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

m2 <- make_comb_mat(hallmarks_m, mode = "intersect")
m2 <- m2[comb_degree(m2) >= 2]
UpSet(m2)
ht = draw(UpSet(m2))
od = column_order(ht)
cs = comb_size(m2)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

m3 <- make_comb_mat(miRNAs, mode = "intersect")
m3 <- m3[comb_degree(m3) >= 2]
UpSet(m3)
ht = draw(UpSet(m3))
od = column_order(ht)
cs = comb_size(m3)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

m4 <- make_comb_mat(truetargets_m, mode = "intersect")
m4 <- m4[comb_degree(m4) >= 2]
UpSet(m4)
ht = draw(UpSet(m4))
od = column_order(ht)
cs = comb_size(m4)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

m5 <- make_comb_mat(dm_m, mode = "intersect")
m5 <- m5[comb_degree(m5) >= 2]
UpSet(m5)
ht = draw(UpSet(m5))
od = column_order(ht)
cs = comb_size(m5)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})


m6 <- make_comb_mat(dmtt_m, mode = "intersect")
m6 <- m6[comb_degree(m6) >= 2]
UpSet(m6)
ht = draw(UpSet(m6))
od = column_order(ht)
cs = comb_size(m6)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})






