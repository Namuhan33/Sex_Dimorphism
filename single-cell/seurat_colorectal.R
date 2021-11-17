setwd("/Users/namuhanz/cancer/single cell/test")
# data preparation ####################################################################
library(GEOquery)
library(R.utils)
library(stringr)
### prepare count matrix of CRC tumor cell and normal mucosal cell
# gunzip("./GSE81861/GSE81861_CRC_tumor_all_cells_COUNT.csv.gz", remove=FALSE)
data.count <- read.csv("./GSE81861/GSE81861_CRC_tumor_all_cells_COUNT.csv", header = T)
genenames <- data.count[,1]
genesymbol <- word(genenames, 2, sep = "_")
genesymbol_unique <- make.names((genesymbol), unique = TRUE)
rownames(data.count) <- genesymbol_unique
data.count <- data.count[,-1]

# gunzip("./GSE81861_CRC_NM_all_cells_COUNT.csv.gz", remove = F)
data.count.nm <- read.csv("./GSE81861_CRC_NM_all_cells_COUNT.csv", header = T)
genenames.nm <- data.count.nm[,1]
genesymbol.nm <- word(genenames.nm, 2, sep = "_")
genesymbol_unique.nm <- make.names((genesymbol.nm), unique = TRUE)
rownames(data.count.nm) <- genesymbol_unique.nm
data.count.nm <- data.count.nm[,-1]
data <- cbind(data.count, data.count.nm)
### prepare sample info matrix
sampleinfo <- read.csv("GSE81861_series_matrix.csv", header = T)
sampleinfo_1 <- sampleinfo[c(37,46,47,48,50),]
rownames(sampleinfo_1) <- c("sample_title", "patient_id", "sample_type", "tissue_type", "patient_gender")
sampleinfo_1 <- sampleinfo_1[,-1]
colnames(sampleinfo_1) <- NULL
patient_id <- as.character(sampleinfo_1["patient_id",])
patient_id <- gsub(".*: ", "", patient_id, perl = TRUE)
tissue_type <- as.character(sampleinfo_1["tissue_type",])
tissue_type <- gsub(".*: ", "", tissue_type, perl = TRUE)
patient_gender <- as.character(sampleinfo_1["patient_gender",])
patient_gender <- gsub(".*: ", "", patient_gender, perl = TRUE)
sampleinfo_2 <- data.frame(as.character(sampleinfo_1["sample_title",]), patient_gender, tissue_type,patient_id)
colnames(sampleinfo_2) <- c("sample_tittle","patient_gender","tissue_type", "patient_id")
# metadata of patients ######################################################################
sampleinfo_3 <- colnames(data)
cellid <- sub("[_].*", "", sampleinfo_3)
celltype <- gsub(".*[_]([^.]+)[_].*", "\\1", sampleinfo_3)
celltype <- sub("[_].*", "", celltype)
sampleinfo_4 <- data.frame(cellid, celltype)
colnames(sampleinfo_2)[1] <- 'cellid'
sample_metadata <- merge(x = sampleinfo_4, y = sampleinfo_2,'cellid')
sample_metadata$sex.tissue <- paste(sample_metadata$patient_gender, sample_metadata$tissue_type, sep = '.')

# sample_metadata <- cbind(sampleinfo_4, sampleinfo_2["sample_tittle" == sampleinfo_4[, "cellid"], ])

# single cell ###############################################################################
library(dplyr)
library(Seurat)
library(patchwork)
### prepare Seurat object
colnames(data) <- substr(colnames(data), 1, 7)
filtered_data <- data[,sample_metadata$cellid]
pbmc <- CreateSeuratObject(counts = filtered_data, project = "SeuratProject", min.cells = 3, min.features = 200)
### set seurat object identity
sample_metadata <- sample_metadata[order(match(sample_metadata[,'cellid'],rownames(pbmc@meta.data))),]
# Idents(object = pbmc) <- sample_metadata$sex.tissue
Idents(object = pbmc) <- sample_metadata$celltype
pbmc@meta.data$cell_type <- sample_metadata$celltype
# pbmc[["cell.ident"]] <- Idents(object = pbmc)
# Idents(object = pbmc) <- sample_metadata$patient_gender
# pbmc[["sex.ident"]] <- Idents(object = pbmc)
# Idents(object = pbmc) <- tissue <- sample_metadata$tissue_type


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT[\\.|\\-]")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
### QC
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 70)
### normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# pbmc <- NormalizeData(pbmc)
### find most variable features/genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(pbmc), 10)
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4
all.genes <- rownames(pbmc)
### linear transformation before PCA (eliminate batch effect as well)
pbmc <- ScaleData(pbmc, features = all.genes)
### PCA summarization (linear dimensional reduction)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
### determine components/dimensionalities based on PCA score
pbmc <- JackStraw(pbmc, num.replicate = 100) ### determine statistical significance of PCA scores
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
Dimplot( obj, group.by = 'cell_type')
# DimPlot(pbmc, reduction = "umap")
# saveRDS(pbmc, file = "./test_tissue.rds")

# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
### min.pct: Pre-filter features that are detected at <50% frequency in selected cells
### min.diff.pct: Pre-filter features whose detection percentages across the two groups are similar (within 0.25) 
# pbmc.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 2, order_by = avg_log2FC)

# subcluster by cell type ##############
newsample_metadata <- subset(sample_metadata, cellid %in% rownames(pbmc@meta.data))
Idents(object = pbmc) <- newsample_metadata$celltype

bcell <- subset(pbmc, idents =  "Bcell")
Idents(object = bcell) <- newsample_metadata[newsample_metadata$celltype=="Bcell",]$sex.tissue

tcell <- subset(pbmc, idents =  "Tcell")
Idents(object = tcell) <- newsample_metadata[newsample_metadata$celltype=="Tcell",]$sex.tissue

macrophage <- subset(pbmc, idents = "Macrophage")
Idents(object = macrophage) <- newsample_metadata[newsample_metadata$celltype=="Macrophage",]$sex.tissue

epithelial <- subset(pbmc, idents = "Epithelial")
Idents(object = epithelial) <- newsample_metadata[newsample_metadata$celltype=="Epithelial",]$sex.tissue

fibroblast <- subset(pbmc, idents = "Fibroblast")
Idents(object = fibroblast) <- newsample_metadata[newsample_metadata$celltype=="Fibroblast",]$sex.tissue

endothelial <- subset(pbmc, idents = "Endothelial")
Idents(object = endothelial) <- newsample_metadata[newsample_metadata$celltype=="Endothelial",]$sex.tissue

mastcell <- subset(pbmc, idents = "MastCell")
Idents(object = mastcell) <- newsample_metadata[newsample_metadata$celltype=="MastCell",]$sex.tissue




female.mastcell <- FindMarkers(mastcell, 
                                 ident.1 = "Female.Colorectal Tumor", 
                                 ident.2 = "Female.Normal Mucosa", 
                                 logfc.threshold = log(1), 
                                 test.use = "wilcox")
female.mastcell_sig <- female.mastcell[female.mastcell$p_val_adj < 0.05,]

male.mastcell <- FindMarkers(mastcell, 
                               ident.1 = "Male.Colorectal Tumor", 
                               ident.2 = "Male.Normal Mucosa",
                               logfc.threshold = log(1), 
                               test.use = "wilcox")
male.mastcell_sig <- male.mastcell[male.mastcell$p_val_adj < 0.05,]


write.csv(female.mastcell_sig, file = "../female_Mastcell.csv")
write.csv(male.mastcell_sig, file = "../male_Mastcell.csv")

# subcluster by sex and tumor or not #########
newsample_metadata <- subset(sample_metadata, cellid %in% rownames(pbmc@meta.data))
Idents(object = pbmc) <- newsample_metadata$sex.tissue

female.pbmc <- FindMarkers(pbmc, 
                               ident.1 = "Female.Colorectal Tumor", 
                               ident.2 = "Female.Normal Mucosa", 
                               logfc.threshold = log(1), 
                               test.use = "wilcox")
female.pbmc_sig <- female.pbmc[female.pbmc$p_val_adj < 0.05,]

male.pbmc <- FindMarkers(pbmc, 
                             ident.1 = "Male.Colorectal Tumor", 
                             ident.2 = "Male.Normal Mucosa",
                             logfc.threshold = log(1), 
                             test.use = "wilcox")
male.pbmc_sig <- male.pbmc[male.pbmc$p_val_adj < 0.05,]

table(newsample_metadata[newsample_metadata$celltype=='Epithelial',]$sex.tissue)

write.csv(female.pbmc_sig, file = "../female_all_de.csv")
write.csv(male.pbmc_sig, file = "../male_all_de.csv")



# Levene's test #######
all_normalized_count <- as.matrix(pbmc@assays$RNA@data)

Idents(object = pbmc) <- newsample_metadata$celltype
epithelial <- subset(pbmc, idents = "Epithelial")
epi_normalized_count <- as.matrix(epithelial@assays$RNA@data)
epi_meta <- newsample_metadata[newsample_metadata$celltype=="Epithelial",]

all_female <- all_normalized_count[,newsample_metadata[newsample_metadata$patient_gender=="Female",]$cellid]
all_male <- all_normalized_count[,newsample_metadata[newsample_metadata$patient_gender=="Male",]$cellid]
epi_female <- epi_normalized_count[,epi_meta[epi_meta$patient_gender=="Female",]$cellid]
epi_male <- epi_normalized_count[,epi_meta[epi_meta$patient_gender=="Male",]$cellid]

all_female_meta <- newsample_metadata[newsample_metadata$patient_gender=="Female",]
all_male_meta <- newsample_metadata[newsample_metadata$patient_gender=="Male",]
epi_female_meta <- epi_meta[epi_meta$patient_gender=="Female",]
epi_male_meta <- epi_meta[epi_meta$patient_gender=="Male",]

all_female_meta <- all_female_meta[,c(1,4)]
all_male_meta <- all_male_meta[,c(1,4)]
epi_female_meta <- epi_female_meta[,c(1,4)]
epi_male_meta <- epi_male_meta[,c(1,4)]

all_female_meta <- all_female_meta[order(all_female_meta$tissue_type, decreasing = T),]
all_male_meta <- all_male_meta[order(all_male_meta$tissue_type, decreasing = T),]
epi_female_meta <- epi_female_meta[order(epi_female_meta$tissue_type, decreasing = T),]
epi_male_meta <- epi_male_meta[order(epi_male_meta$tissue_type, decreasing = T),]


all_female <- all_female[,all_female_meta$cellid]
all_male <- all_male[,all_male_meta$cellid]
epi_female <- epi_female[,epi_female_meta$cellid]
epi_male <- epi_male[,epi_male_meta$cellid]

library(car)
### all female #########
CountsforLevene_af <- data.frame()
Tissue_list <- as.factor(c(rep(1,137*nrow(all_female)), rep(2,186*nrow(all_female))))
Gene_list <- data.frame(c(rep(rownames(all_female),137), rep(rownames(all_female),186)))
Counts_list <- data.frame(x=unlist(as.data.frame(all_female)))
CountsforLevene_af <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene_af) <- c("Tissue", "Gene", "Counts")
n <- nrow(all_female)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(all_female)[i]
  CountsforLevene_gene <- CountsforLevene_af[CountsforLevene_af$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "../DV/DV_all_all_female.csv")
write.csv(Levene_results_sig, file = "../DV/DV_sig_all_female.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- all_female[DVgene,]
nor <- character()
tum <- character()
for(a in 1:length(DVgene)){
  normal <- as.numeric(DVcount[a, c(1:137)])
  tumor <- as.numeric(DVcount[a,c(137+1:186)])
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
write.csv(DV_MAD, file = "../DV/DV_MAD_all_female.csv")
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "../DV/result_all_female.csv")

### all male #########
CountsforLevene_am <- data.frame()
Tissue_list <- as.factor(c(rep(1,59*nrow(all_male)), rep(2,120*nrow(all_male))))
Gene_list <- data.frame(c(rep(rownames(all_male),59), rep(rownames(all_male),120)))
Counts_list <- data.frame(x=unlist(as.data.frame(all_male)))
CountsforLevene_am <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene_am) <- c("Tissue", "Gene", "Counts")
n <- nrow(all_male)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(all_male)[i]
  CountsforLevene_gene <- CountsforLevene_am[CountsforLevene_am$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "../DV/DV_all_all_male.csv")
write.csv(Levene_results_sig, file = "../DV/DV_sig_all_male.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- all_male[DVgene,]
nor <- character()
tum <- character()
for(a in 1:length(DVgene)){
  normal <- as.numeric(DVcount[a, c(1:59)])
  tumor <- as.numeric(DVcount[a,c(59+1:120)])
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
write.csv(DV_MAD, file = "../DV/DV_MAD_all_male.csv")
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "../DV/result_all_male.csv")

### epi female #########
CountsforLevene_ef <- data.frame()
Tissue_list <- as.factor(c(rep(1,101*nrow(epi_female)), rep(2,126*nrow(epi_female))))
Gene_list <- data.frame(c(rep(rownames(epi_female),101), rep(rownames(epi_female),126)))
Counts_list <- data.frame(x=unlist(as.data.frame(epi_female)))
CountsforLevene_ef <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene_ef) <- c("Tissue", "Gene", "Counts")
n <- nrow(epi_female)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(epi_female)[i]
  CountsforLevene_gene <- CountsforLevene_ef[CountsforLevene_ef$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "../DV/DV_all_epi_female.csv")
write.csv(Levene_results_sig, file = "../DV/DV_sig_epi_female.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- epi_female[DVgene,]
nor <- character()
tum <- character()
for(a in 1:length(DVgene)){
  normal <- as.numeric(DVcount[a, c(1:101)])
  tumor <- as.numeric(DVcount[a,c(101+1:126)])
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
write.csv(DV_MAD, file = "../DV/DV_MAD_epi_female.csv")
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "../DV/result_epi_female.csv")


### epi male #########
CountsforLevene_em <- data.frame()
Tissue_list <- as.factor(c(rep(1,40*nrow(epi_male)), rep(2,78*nrow(epi_male))))
Gene_list <- data.frame(c(rep(rownames(epi_male),40), rep(rownames(epi_male),78)))
Counts_list <- data.frame(x=unlist(as.data.frame(epi_male)))
CountsforLevene_em <- cbind(Tissue_list,Gene_list,Counts_list)
colnames(CountsforLevene_em) <- c("Tissue", "Gene", "Counts")
n <- nrow(epi_male)
Levene_results <- data.frame( Gene_Name = character(n), 
                              P_Value = numeric(n), 
                              F_Value = numeric(n), 
                              stringsAsFactors = FALSE)
for(i in 1:n){
  gene <- rownames(epi_male)[i]
  CountsforLevene_gene <- CountsforLevene_em[CountsforLevene_em$Gene==gene,]
  res.levenetest <- leveneTest(Counts ~ Tissue, data = CountsforLevene_gene)
  Levene_results$Gene_Name[i] <- c(gene)
  Levene_results$P_Value[i] <- c(res.levenetest$'Pr(>F)'[1])
  Levene_results$F_Value[i] <- c(res.levenetest$'F value'[1])
}
adjP <- p.adjust(Levene_results$P_Value, method = "BH", n = length(Levene_results$P_Value))
Levene_results <- cbind(Levene_results, adjP)
Levene_results_sig <- subset(Levene_results, adjP < 0.05)
write.csv(Levene_results, file = "../DV/DV_all_epi_male.csv")
write.csv(Levene_results_sig, file = "../DV/DV_sig_epi_male.csv")
DVgene <- Levene_results_sig$Gene_Name
DVcount <- epi_male[DVgene,]
nor <- character()
tum <- character()
for(a in 1:length(DVgene)){
  normal <- as.numeric(DVcount[a, c(1:40)])
  tumor <- as.numeric(DVcount[a,c(40+1:78)])
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
write.csv(DV_MAD, file = "../DV/DV_MAD_epi_male.csv")
result <- aggregate(data.frame(count = type), list(value = type), length)
write.csv(result, file = "../DV/result_epi_male.csv")

