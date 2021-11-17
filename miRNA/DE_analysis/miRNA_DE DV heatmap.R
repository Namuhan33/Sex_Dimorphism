### input counts data ###############################################################
setwd("/Users/namuhanz/cancer/miRNA/DE DV heatmap density/file")
countsdata <- read.csv("./Counts_homo.csv", header = FALSE)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)

num_n <- 19 ### input the number of normal sample
num_t <- 21 ### input the number of tumor sample
### filtering and normalization ########################################################
library(edgeR)
group <- factor(c(rep("normal",num_n),rep("tumor",num_t))) 
y <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

normalizedcounts <- y$counts

design <- model.matrix(~group)
y <- estimateDisp(y)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
fdr <- summary(decideTests(qlf))
qlf <- topTags(qlf, n=100000)
qlf <- as.data.frame(qlf)
qlf <- cbind(rownames(qlf), qlf)
colnames(qlf) <- c("gene_id", "log2FoldChange", "log2CPM", "LR","PValue", "FDR")
write.table(qlf, "./DEmiRNA_all.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlf_FDR <- qlf[which(qlf$FDR <= 0.05), ]
qlf_FDR[which(qlf_FDR$log2FoldChange > 0), "up-down"] <- "Up"
qlf_FDR[which(qlf_FDR$log2FoldChange < 0), "up-down"] <- "Down"
write.table(qlf_FDR, "./DEmiRNA_FDR.xls", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")
qlfsig <- qlf[which(qlf$FDR <= 0.05 & abs(qlf$log2FoldChange) >= 1) , ]
qlfsig[which(qlfsig$log2FoldChange > 0), "up-down"] <- "Up"
qlfsig[which(qlfsig$log2FoldChange < 0), "up-down"] <- "Down"
write.table(qlfsig, "./DEmiRNA_FDR&FC.xls", sep ="\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
numofup <- length(which(qlfsig$'up-down' == "Up"))
numofdown <- length(which(qlfsig$'up-down' == "Down"))
### heatmap ######################################################################
library(RColorBrewer)
library(pheatmap)
library(viridis)

nr <- nrow(qlfsig)
nc <- ncol(normalizedcounts)
degcounts <- data.frame(matrix(NA,
                               nrow = nr,
                               ncol = nc))
colnames(degcounts) <- colnames(normalizedcounts)
rownames(degcounts) <- rownames(qlfsig)

for(gene in qlfsig$gene_id){
  degcounts[gene,] <- normalizedcounts[gene,]
}

degcounts_log <- cpm(degcounts, log = T, prior.count = 1)
Z_degcounts_log <- t(scale(t(degcounts_log)))

col_groups <- group

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(degcounts)
mat_colors <- list(group = brewer.pal(3, "Set2"))
names(mat_colors$group) <- unique(col_groups)
pheatmap(
  mat               = Z_degcounts_log,
  scale             = "row",
  color             = inferno(20),
  border_color      = NA,
  show_colnames     = F,
  show_rownames     = T,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "DE miRNAs Heatmap"
)
### density plot ################################################################################
# library("openanalytics/ganalyse")
# x <- DGEList(counts = countdata, group = group)
# density_plot(x, 
#              interactive = FALSE, 
#              title = "Density Plot of miRNA Raw Counts", 
#              groups = group, 
#              facet_cols = 1L,
#              log = FALSE, 
#              # filename, 
#              height = 8L, 
#              width = 12L)







### DV ###########################################################################
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_n*nrow(normalizedcounts)), rep(2,num_t*nrow(normalizedcounts))))
Gene_list <- data.frame(c(rep(rownames(normalizedcounts),num_n), rep(rownames(normalizedcounts),num_t)))
Counts_list <- data.frame(x=unlist(as.data.frame(normalizedcounts)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")

library(car)
n <- nrow(normalizedcounts)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)

for(i in 1:n){
  gene <- rownames(normalizedcounts)[i]
  CountsforLevene_gene <- CountsforLevene[CountsforLevene$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}

adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP <= 0.05)

write.csv(Levene_results, file = "./Var_all.csv")
write.csv(Levene_results_sig, file = "./Var_sig.csv")

summary(keep)
print(fdr)
print(c(numofup, numofdown))

print(nrow(Levene_results_sig))














