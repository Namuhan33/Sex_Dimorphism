setwd("/Users/namuhanz/cancer/miRNA/DE DV heatmap density/file")
countsdata <- read.csv("./Counts_homo.csv", header = FALSE)
colnames(countsdata) <- countsdata[1,]
rownames(countsdata) <- countsdata[,1]
countsdata <- countsdata[-1,]
countsdata <- countsdata[,-1]
countdata <- apply(as.matrix(countsdata), 2, as.numeric)
rownames(countdata) <- rownames(countsdata)
colnames(countdata) <- substr(colnames(countdata),start=9,stop=16)

qlfsig <- read.table("DEmiRNA_FDR&FC.xls", header = T)

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
counts.cpm <- cpm(y, log = F)
counts.log <- log(counts.cpm+1)
Z_counts_log <- t(scale(t(counts.log)))
normalizedcounts <- Z_counts_log


library(RColorBrewer)
library(pheatmap)
library(viridis)

nr <- nrow(qlfsig)
nc <- ncol(normalizedcounts)
degcounts <- data.matrix(matrix(NA,
                               nrow = nr,
                               ncol = nc))
colnames(degcounts) <- colnames(normalizedcounts)
rownames(degcounts) <- qlfsig[,1]

for(gene in qlfsig$gene_id){
  degcounts[gene,] <- normalizedcounts[gene,]
}



col_groups <- group

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(degcounts)
mat_colors <- list(group = RColorBrewer::brewer.pal(2, "Set2"))
names(mat_colors$group) <- unique(col_groups)
pheatmap(
  mat               = degcounts,
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

### DV ###########################################################################
CountsforLevene <- data.frame()
Tissue_list <- as.factor(c(rep(1,num_n*nrow(counts.log)), rep(2,num_t*nrow(counts.log))))
Gene_list <- data.frame(c(rep(rownames(counts.log),num_n), rep(rownames(counts.log),num_t)))
Counts_list <- data.frame(x=unlist(as.data.frame(counts.log)))
CountsforLevene <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene) <- c("Tissue", "Gene", "Counts")

library(car)
n <- nrow(counts.log)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)

for(i in 1:n){
  gene <- rownames(counts.log)[i]
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



print(nrow(Levene_results_sig))


