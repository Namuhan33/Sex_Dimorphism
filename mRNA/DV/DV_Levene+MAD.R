###########################################################################################
### input counts matrix
setwd("/Users/namuhanz/cancer/transcripts/expression variability/file")
countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]

countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)

num_nor <- 58 ### input the number of normal samples
num_tum <- 58 ### inpout the number of tumor samples
###########################################################################################
### filter out lowly expressed genes and normalize the library size of the counts
library(edgeR)
group <- factor(c(rep("normal",num_nor),rep("tumor",num_tum))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
###########################################################################################
### create data.frame of counts for Levene-test
counts.cpm <- cpm(y,log=F)
counts.log <- log(counts.cpm + 1)
normalized_countdata <- counts.log
CountsforLevene <- data.frame()

Tissue_list <- as.factor(c(rep(1,num_nor*nrow(normalized_countdata)), rep(2,num_tum*nrow(normalized_countdata))))
Gene_list <- data.frame(c(rep(rownames(normalized_countdata),num_nor), rep(rownames(normalized_countdata),num_tum)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalized_countdata)))

CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")

###########################################################################################
### Levene's test on genes one-by-one
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
# adjP_1 <- p.adjust(Levene_results$P_Value, method = "fdr", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)


# Variability_sig <- subset(Levene_results_sig, Gene_Name != 0 & P_Value != 0 & F_Value != 0)
write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")

# MAD ##################
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





