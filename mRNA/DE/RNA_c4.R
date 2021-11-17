### read in counts matrix ##########################################################################################
setwd("/Users/namuhanz/cancer/transcripts/disease/file")
library(edgeR)
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
num_mn <- 10 ### input the number of male normal sample
num_mt <- 12 ### input the number of male tumor sample
num_fn <- 9 ### input the number of female normal sample
num_ft <- 9 ### input the number of female tumor sample

### setup sample description data frame #############################################################################################
sampledescription <- data.frame()
SampleID <- data.frame(colnames(countdata))
Sex <- data.frame(c(rep("male", num_mn), 
                    rep("male", num_mt), 
                    rep("female", num_fn), 
                    rep("female", num_ft)))
Type <- data.frame(c(rep("normal", num_mn), 
                     rep("tumor", num_mt), 
                     rep("normal", num_fn), 
                     rep("tumor", num_ft)))
sampledescription <- cbind(SampleID, Sex, Type)
colnames(sampledescription) <- c("SampleID", "Sex", "Type")
### DE analysis #############################################################################################
### create design matrix
sex <- factor(c(rep("male", num_mn), 
                 rep("male", num_mt), 
                 rep("female", num_fn), 
                 rep("female", num_ft)))
type <- factor(c(rep("normal", num_mn), 
                  rep("tumor", num_mt), 
                  rep("normal", num_fn), 
                  rep("tumor", num_ft)))
design <- model.matrix(~type+sex)

### create DEG list
group <- factor(paste0(sampledescription$Sex, ".", sampledescription$Type))
y <- DGEList(countdata, group = group)
rownames(design) <- colnames(y)
### filtering and normalization
keep <- filterByExpr(y,design)
summary(keep)
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
### estimating the dispersion
y <- estimateDisp(y, design, robust = T)
### differential expression
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = ncol(fit$design))
fdr <- summary(decideTests(qlf))
### output
qlf <- topTags(qlf, n=100000)
qlf <- as.data.frame(qlf)
qlf <- cbind(rownames(qlf), qlf)

colnames(qlf) <- c("gene_id", "log2FoldChange", "log2CPM", "LR","PValue", "FDR")
write.table(qlf, "./DEG_all.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlf_FDR <- qlf[which(qlf$FDR < 0.05), ]
qlf_FDR[which(qlf_FDR$log2FoldChange > 0), "up-down"] <- "Up"
qlf_FDR[which(qlf_FDR$log2FoldChange < 0), "up-down"] <- "Down"
write.table(qlf_FDR, "./DEG_FDR.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlfsig <- qlf[which(qlf$FDR < 0.05 & abs(qlf$log2FoldChange) >= 1) , ]
qlfsig[which(qlfsig$log2FoldChange > 0), "up-down"] <- "Up"
qlfsig[which(qlfsig$log2FoldChange < 0), "up-down"] <- "Down"
numofup <- length(which(qlfsig$'up-down' == "Up"))
numofdown <- length(which(qlfsig$'up-down' == "Down"))
write.table(qlfsig, "./DEG_FDR&FC.xls", sep ="\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")

print(summary(keep))
print(fdr)
print(c(numofup, numofdown))






