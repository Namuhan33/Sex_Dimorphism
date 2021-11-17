setwd("/Users/namuhanz/cancer/methylation/exact sample names/file")
library(rjson)
###get a df matrix: replace file names with Barcodes
metadata_json_File <- fromJSON(file = "./metadata.cart.2021-09-17.json")
json_File_Info <-data.frame(filesName = c(), TCGA_Barcode = c())
for(i in 1:length(metadata_json_File)){
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <- metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info, data.frame(filesName = file_name, TCGA_Barcode = TCGA_Barcode))
}
rownames(json_File_Info) <- json_File_Info[,1]

library(TCGAbiolinks)
query <- GDCquery(project = c("TCGA-THCA","TCGA-DLBC"),
                  data.category = "DNA Methylation",
                  data.type = "Methylation Beta Value",
                  workflow.type = "Liftover",
                  platform = "Illumina Human Methylation 450") #change every analysis
SamplesDown <- getResults(query, cols = c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = SamplesDown, typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = SamplesDown, typesample = "NT")

library(dplyr)
a <- intersect(dataSmNT, json_File_Info$TCGA_Barcode)
b <- intersect(dataSmTP, json_File_Info$TCGA_Barcode)

df <- data.frame(matrix(NA,
                        ncol = length(a)+length(b),
                        nrow = 2))
df[1,] <- c(a,b)
colnames(df) <- df[1,]
df[2,] <- get_IDs(df)[,3]
df <- df[-1,]




df[nrow(df)+1,] <- 0

H <- new.env()
H_1 <- new.env()
H_2 <- new.env()

for (i in 1:ncol(df)) {
  H[[df[1,i]]] <-  0
  H_1[[df[1,i]]] <-  0
  H_2[[df[1,i]]] <-  0
}

for (i in 1:ncol(df)) {
  H[[df[1,i]]] <- H[[df[1,i]]] + 1
}

for (i in 1:length(a)) { 
  H_1[[df[1,i]]] <- H_1[[df[1,i]]] + 1
}

for (i in length(a)+1 :ncol(df)) {
  H_2[[df[1,i]]] <- H_2[[df[1,i]]] + 1
}

for (i in 1:ncol(df)) {
  if (H[[df[1,i]]] > 1 & H_1[[df[1,i]]] >= 1 & H_2[[df[1,i]]] >= 1) {
    # if (H[[df[1,i]]] > 1) {
    df[nrow(df),i] <- 1
  }
}

df <- df[,df[nrow(df),] == 1]

df <- df[3:nrow(df)-1,]

num_of_normal <- length(intersect(dataSmNT, colnames(df)))
num_of_tumor <- length(intersect(dataSmTP, colnames(df)))

samplenames <- colnames(df)
write.csv(samplenames, file = "./samplenames.csv")





