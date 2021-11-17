setwd("/Users/namuhanz/cancer/methylation/exact sample names/file")
RNAcounts <- read.csv("./Counts_homo.csv", header = F)
num_nor <- 17
RNAcounts <- RNAcounts[,-1]
colnames(RNAcounts) <- RNAcounts[1,]
RNAsamples <- colnames(RNAcounts)[c(1:num_nor)]
samplenames <- character()
for(i in 1:num_nor){
  samplenames[i] <- substr(RNAsamples[i],1,12)
}
write.csv(samplenames, file = "./samplenames.csv", row.names = F)