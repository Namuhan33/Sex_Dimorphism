### the interaction of sex within type
setwd("/Users/namuhanz/cancer/transcripts/disease/file")
library(edgeR)
num_mn <- 10 ### input the number of male normal sample
num_mt <- 12 ### input the number of male tumor sample
num_fn <- 9 ### input the number of female normal sample
num_ft <- 9 ### input the number of female tumor sample
#data prep ############
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

#description matrix ######
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


target_1 <- target
target_1$Sex <- relevel(target_1$Sex, ref = "male")

#DGE list create ######
group <- factor(paste0(sampledescription$Sex, ".", sampledescription$Type))
y <- DGEList(countdata, group = group)
y_1 <- DGEList(countdata, group = group)

#design matrix #######
design <- model.matrix(~Sex+Type+Sex:Type, data = target)
design_1 <- model.matrix(~Sex+Type+Sex:Type, data=target_1)
### equal to ~sex*type, data = target
#filtering and normalization #####
keep <- filterByExpr(y,design)
keep_1 <- filterByExpr(y_1, design_1)

summary(keep)

y <- y[keep, ,keep.lib.sizes=FALSE]
y_1 <- y_1[keep_1, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y_1 <- calcNormFactors(y_1)
y <- estimateDisp(y, design, robust = T)
y_1 <- estimateDisp(y_1,design_1, robust = T)
fit <- glmQLFit(y, design)
fit_1 <- glmQLFit(y_1, design_1)

# contrasts and compare #########
qlf <- glmQLFTest(fit, coef = 4)
summary(decideTests(qlf))
fdr <- summary(decideTests(qlf))
qlf_1 <- glmQLFTest(fit_1,coef = 4)
summary(decideTests(qlf_1))
fdr_1 <- summary(decideTests(qlf_1))

# tidy up results #########
qlf <- topTags(qlf, n=100000)
qlf <- as.data.frame(qlf)
qlf <- cbind(rownames(qlf), qlf)
colnames(qlf) <- c("gene_id", "log2FoldChange", "log2CPM", "F","PValue", "FDR")
write.table(qlf, "./DEG_all_m.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlf_FDR <- qlf[which(qlf$FDR <= 0.05), ]
qlf_FDR[which(qlf_FDR$log2FoldChange > 0), "up-down"] <- "Up"
qlf_FDR[which(qlf_FDR$log2FoldChange < 0), "up-down"] <- "Down"
write.table(qlf_FDR, "./DEG_FDR_m.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlfsig <- qlf[which(qlf$FDR < 0.05 & abs(qlf$log2FoldChange) >= 1) , ]
qlfsig[which(qlfsig$log2FoldChange > 0), "up-down"] <- "Up"
qlfsig[which(qlfsig$log2FoldChange < 0), "up-down"] <- "Down"
numofup <- length(which(qlfsig$'up-down' == "Up"))
numofdown <- length(which(qlfsig$'up-down' == "Down"))
write.table(qlfsig, "./DEG_FDR&FC_m.xls", sep ="\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")


qlf_1 <- topTags(qlf_1, n=100000)
qlf_1 <- as.data.frame(qlf_1)
qlf_1 <- cbind(rownames(qlf_1), qlf_1)
colnames(qlf_1) <- c("gene_id", "log2FoldChange", "log2CPM", "F","PValue", "FDR")
write.table(qlf_1, "./DEG_all_f.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlf_1_FDR <- qlf_1[which(qlf_1$FDR < 0.05), ]
qlf_1_FDR[which(qlf_1_FDR$log2FoldChange > 0), "up-down"] <- "Up"
qlf_1_FDR[which(qlf_1_FDR$log2FoldChange < 0), "up-down"] <- "Down"
write.table(qlf_1_FDR, "./DEG_FDR_f.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlf_1sig <- qlf_1[which(qlf_1$FDR < 0.05 & abs(qlf_1$log2FoldChange) >= 1) , ]
qlf_1sig[which(qlf_1sig$log2FoldChange > 0), "up-down"] <- "Up"
qlf_1sig[which(qlf_1sig$log2FoldChange < 0), "up-down"] <- "Down"
numofup_1 <- length(which(qlf_1sig$'up-down' == "Up"))
numofdown_1 <- length(which(qlf_1sig$'up-down' == "Down"))
write.table(qlf_1sig, "./DEG_FDR&FC_f.xls", sep ="\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")


# take down the results #######
print(summary(keep))
print(fdr)
print(fdr_1)
print(c(numofup, numofdown))
print(c(numofup_1, numofdown_1))




