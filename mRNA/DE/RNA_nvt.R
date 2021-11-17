setwd("/Users/namuhanz/cancer/transcripts/disease/file")
library(edgeR)

countsdata <- read.csv("./Counts_homo.csv", header = F)
colnames(countsdata) <- countsdata[1,]
countsdata <- countsdata[-1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)

num_n <- 17
num_t <- 17

group <- factor(c(rep("normal",num_n),rep("tumor",num_t))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~group)
y <- estimateDisp(y)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
fdr <- summary(decideTests(qlf))
qlf <- topTags(qlf, n=100000)
qlf <- as.data.frame(qlf)
qlf <- cbind(rownames(qlf), qlf)
colnames(qlf) <- c("gene_id", "log2FoldChange", "log2CPM", "LR","PValue", "FDR")
write.table(qlf, "./DE_all.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlf_FDR <- qlf[which(qlf$FDR < 0.05), ]
qlf_FDR[which(qlf_FDR$log2FoldChange > 0), "up-down"] <- "Up"
qlf_FDR[which(qlf_FDR$log2FoldChange < 0), "up-down"] <- "Down"
write.table(qlf_FDR, "./DE_FDR.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlfsig <- qlf[which(qlf$FDR < 0.05 & abs(qlf$log2FoldChange) >= 1) , ]
qlfsig[which(qlfsig$log2FoldChange > 0), "up-down"] <- "Up"
qlfsig[which(qlfsig$log2FoldChange < 0), "up-down"] <- "Down"
write.table(qlfsig, "./DE_FDR&FC.xls", sep ="\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
numofup <- length(which(qlfsig$'up-down' == "Up"))
numofdown <- length(which(qlfsig$'up-down' == "Down"))

summary(keep)
print(fdr)
print(c(numofup, numofdown))
