### move all counts txt files to one file
setwd("/Users/namuhanz/cancer/miRNA/raw file processing")
options(stringsAsFactors = F )
dir.create("SampleFiles")
filepath <- dir (path = "./miRNA Expression Quantification", full.names = TRUE)
for(wd in filepath){
  files <- dir(path = wd, pattern = "txt$")
  fromfilepath <- paste(wd, "/", files, sep = "")
  tofilepath <- paste("./SampleFiles/", files, sep = "" )
  file.copy(fromfilepath, tofilepath)
}

setwd("./SampleFiles")
file.remove("annotations.txt")
file.remove("MANIFEST.txt")



###replace file names with Barcodes
setwd("/Users/namuhanz/cancer/miRNA/raw file processing")
library(rjson)
metadata_json_File <- fromJSON(file = "./metadata.cart.2021-09-02.json")
json_File_Info <-data.frame(filesName = c(), TCGA_Barcode = c())
for(i in 1:length(metadata_json_File)){
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <- metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info, data.frame(filesName = file_name, TCGA_Barcode = TCGA_Barcode))
}
rownames(json_File_Info) <- json_File_Info[,1]
#write.csv(json_File_Info, file = "./json_File_Info.csv")
filesName_To_TCGA_BarcodeFile <- json_File_Info[-1]
setwd("./SampleFiles")
countsFileNames <- dir(pattern = "txt$")

### replace the column names with barcodes of each file
allsampleRawCounts <- data.frame()
for(txtFile in countsFileNames){
  SampleCounts <- read.table(txtFile, header = T)
  rownames(SampleCounts) <- SampleCounts[,1]
  SampleCounts <- SampleCounts[3]
  colnames(SampleCounts) <- filesName_To_TCGA_BarcodeFile[c(txtFile),]
  if (dim(allsampleRawCounts)[1]== 0){
    allsampleRawCounts <- SampleCounts
  }
  else
  {
    allsampleRawCounts <- cbind(allsampleRawCounts, SampleCounts)
  }
}
# write.csv(allsampleRawCounts, file = "../RawCounts.csv")

### distinguish the samples by tissue type
library(TCGAbiolinks)

query <- GDCquery(project = c("TCGA-THCA","TCGA-DLBC"),
                  data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification",
                  workflow.type = "BCGSC miRNA Profiling") ### change every analysis

SamplesDown <- getResults(query, cols = c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = SamplesDown, typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = SamplesDown, typesample = "NT")

library(dplyr)
a <- intersect(dataSmNT, colnames(allsampleRawCounts))
b <- intersect(dataSmTP, colnames(allsampleRawCounts))

Counts <- data.frame(c(allsampleRawCounts[,a],allsampleRawCounts[,b]))
rownames(Counts) <- rownames(allsampleRawCounts)
colnames(Counts) <- c(a, b)

print(length(a))
print(length(b))

### select out the paired samples of normal and tumor type
Counts <- rbind(get_IDs(Counts)[,3], Counts)
Counts[nrow(Counts)+1,] <- 0

H <- new.env()
H_1 <- new.env()
H_2 <- new.env()

for (i in 1:ncol(Counts)) {
  H[[Counts[1,i]]] <-  0
  H_1[[Counts[1,i]]] <-  0
  H_2[[Counts[1,i]]] <-  0
}

for (i in 1:ncol(Counts)) {
  H[[Counts[1,i]]] <- H[[Counts[1,i]]] + 1
}

for (i in 1:length(a)) { 
  H_1[[Counts[1,i]]] <- H_1[[Counts[1,i]]] + 1
}

for (i in length(a)+1 :ncol(Counts)) {
  H_2[[Counts[1,i]]] <- H_2[[Counts[1,i]]] + 1
}

for (i in 1:ncol(Counts)) {
  if (H[[Counts[1,i]]] > 1 & H_1[[Counts[1,i]]] >= 1 & H_2[[Counts[1,i]]] >= 1) {
    # if (H[[Counts[1,i]]] > 1) {
    Counts[nrow(Counts),i] <- 1
  }
}

Counts <- Counts[,Counts[nrow(Counts),] == 1]

Counts <- Counts[3:nrow(Counts)-1,]

num_of_normal <- length(intersect(dataSmNT, colnames(Counts)))
num_of_tumor <- length(intersect(dataSmTP, colnames(Counts)))

write.csv(Counts, file = "../Counts_homo.csv")







