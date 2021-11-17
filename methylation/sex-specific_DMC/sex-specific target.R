setwd("/Users/namuhanz/cancer/methylation/DM comparison 2vs3/file")
cgidconvert <- read.csv("../cgconvert.csv", header = T)
cgidconvert <- cgidconvert[,-1]
femaledm_df <- read.csv("./female_specific_DM.csv", header = T)
maledm_df <- read.csv("./male_specific_DM.csv", header = T)
femaledm_hyper <- femaledm_df[femaledm_df$sig=="hyper",]$cgid
femaledm_hypo <- femaledm_df[femaledm_df$sig=="hypo",]$cgid
maledm_hyper <- maledm_df[maledm_df$sig=="hyper",]$cgid
maledm_hypo <- maledm_df[maledm_df$sig=="hypo",]$cgid

dmlist <- list(femaledm_hyper,femaledm_hypo,maledm_hyper,maledm_hypo)
targetlist <- list(femaledm_hyper_target=character(),
                   femaledm_hypo_target=character(),
                   maledm_hyper_target=character(),
                   maledm_hypo_target=character())

for (i in 1:length(dmlist)) {
  for (j in 1:length(dmlist[[i]])) {
    cg <- dmlist[[i]][[j]]
    targetlist[[i]][[j]] <- cgidconvert[cgidconvert$cg_id==cg,]$gene_symbols
  }
}

for (k in 1:length(targetlist)) {
  targetlist[[k]] <- targetlist[[k]][targetlist[[k]]!="."]
  targetlist[[k]] <- unlist(strsplit(targetlist[[k]], split = ", "))
  targetlist[[k]] <- unique(targetlist[[k]])
}

for (i in 1:length(targetlist)) {
  path <- paste(names(targetlist[i]), ".csv", sep = "")
  write.csv(targetlist[[i]], file = path)
}





