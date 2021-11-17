setwd("/Users/namuhanz/cancer/transcripts/expression variability/file")
library(edgeR)
num_nor <- 19 ### input the number of normal samples
num_tum <- 21 ### inpout the number of tumor samples
### data read #######
DVlist <- read.csv("./Var_sig.csv", header = T)
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
### filtering and normalization ########
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
normalized_countdata <- y$counts
### select out DV genes ######
DVgene <- DVlist$Gene_Name
DVcount <- normalized_countdata[DVgene,]
### MAD calculate #######
nor <- character()
tum <- character()
for(a in 1:length(DVgene)){
  normal <- as.numeric(DVcount[a, c(1:num_nor)])
  tumor <- as.numeric(DVcount[a,c(num_nor+1:num_tum)])
  nor[a] <- mad(normal)
  tum[a] <- mad(tumor)
}

DV_MAD <- cbind(DVgene,nor,tum)

type <- character()
for(b in 1:length(DVgene)){
  if(as.numeric(DV_MAD[b,2]) > as.numeric(DV_MAD[b,3])){
    type[b] <- "N"
  }
  else{
    type[b] <- "T"
  }
}

DV_MAD <- cbind(DV_MAD,type)
write.csv(DV_MAD, file = "DV_MAD.csv")
aggregate(data.frame(count = type), list(value = type), length)
