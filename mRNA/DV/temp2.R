# bladder all #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/bladder/all")
num_nor <- 19 ### input the number of normal samples
num_tum <- 21 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# bladder female #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/bladder/female")
num_nor <- 9 ### input the number of normal samples
num_tum <- 9 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())





# bladder male #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/bladder/male")
num_nor <- 10 ### input the number of normal samples
num_tum <- 12 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())





# colon all #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/colon/all")
num_nor <- 42 ### input the number of normal samples
num_tum <- 47 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# colon female #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/colon/female")
num_nor <- 22 ### input the number of normal samples
num_tum <- 25 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# colon male #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/colon/male")
num_nor <- 20 ### input the number of normal samples
num_tum <- 22 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# kidney all #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/kidney/all")
num_nor <- 126 ### input the number of normal samples
num_tum <- 126 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# kidney female #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/kidney/female")
num_nor <- 41 ### input the number of normal samples
num_tum <- 41 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# kidney male #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/kidney/male")
num_nor <- 85 ### input the number of normal samples
num_tum <- 85 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# liver all #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/liver/all")
num_nor <- 58 ### input the number of normal samples
num_tum <- 58 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# liver female #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/liver/female")
num_nor <- 25 ### input the number of normal samples
num_tum <- 25 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# liver male #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/liver/male")
num_nor <- 33 ### input the number of normal samples
num_tum <- 33 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# stomach all #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/stomach/all")
num_nor <- 27 ### input the number of normal samples
num_tum <- 27 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# stomach female #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/stomach/female")
num_nor <- 7 ### input the number of normal samples
num_tum <- 7 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# stomach male #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/stomach/male")
num_nor <- 20 ### input the number of normal samples
num_tum <- 20 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# thyroid all #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/thyroid/all")
num_nor <- 58 ### input the number of normal samples
num_tum <- 58 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# thyroid female #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/thyroid/female")
num_nor <- 41 ### input the number of normal samples
num_tum <- 41 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())

# thyroid male #######
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file/thyroid/male")
num_nor <- 17 ### input the number of normal samples
num_tum <- 17 ### inpout the number of tumor samples
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")
library(car)
n <- nrow(normalized_countdata)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(normalized_countdata)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- normalized_countdata[DVgene,]
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
summary(keep)
print(nrow(Levene_results_sig))
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "./result.csv")

rm(list=ls())
