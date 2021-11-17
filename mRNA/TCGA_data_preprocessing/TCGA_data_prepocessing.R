######################################################################################################
###set working directory in the file containing all data files
###rename the file name of the unzipped directory as Gene_Expression_Quantification
setwd("/Users/namuhanz/cancer/transcripts/disease")
options(stringsAsFactors = F )
dir.create("SampleFiles")
filepath <- dir (path = "./Gene_Expression_Quantification", full.names = TRUE)
for(wd in filepath){
  files <- dir(path = wd, pattern = "gz$")
  fromfilepath <- paste(wd, "/", files, sep = "")
  tofilepath <- paste("./SampleFiles/", files, sep = "" )
  file.copy(fromfilepath, tofilepath)
}

setwd("./SampleFiles")

library(R.utils)
countsFiles <- dir(path = "./", pattern = "gz$")
sapply(countsFiles, gunzip) 
###unzip all files and delete original files
###process pattern = json file
###set working directory back to the files containing all data

#######################################################################################################
setwd("/Users/namuhanz/cancer/transcripts/disease") 

library(rjson)
###get a counts matrix: replace file names with Barcodes
metadata_json_File <- fromJSON(file = "./metadata.cart.2021-08-20.json")
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
countsFileNames <- dir(pattern = "counts$")

#######################################################################################################
allsampleRawCounts <- data.frame()
for(txtFile in countsFileNames){
  SampleCounts <- read.table(txtFile, header = FALSE)
  rownames(SampleCounts) <- SampleCounts[,1]
  SampleCounts <- SampleCounts[-1]
  colnames(SampleCounts) <- filesName_To_TCGA_BarcodeFile[paste(txtFile,".gz",sep = ""), ]
  if (dim(allsampleRawCounts)[1]== 0){
    allsampleRawCounts <- SampleCounts
  }
  else
  {
    allsampleRawCounts <- cbind(allsampleRawCounts, SampleCounts)
  }
}
# write.csv(allsampleRawCounts,file = "../allSampleRawCounts.csv")
ensembl_id <- substr(row.names(allsampleRawCounts), 1, 15)
rownames(allsampleRawCounts) <- ensembl_id
# write.csv(allsampleRawCounts, file = "../RawCounts.csv")


RawCounts <- allsampleRawCounts
ensembl_id <- data.frame(ensembl_id = row.names(RawCounts))
rownames(ensembl_id) <- ensembl_id[,1]
RawCounts <- cbind(ensembl_id, RawCounts)

#######################################################################################################
###use gtf file to gain the gene names by ensembl id
get_map = function(input) {
  if(is.character(input)) {
    if(!file.exists(input)) stop("Bad input file.")
    message("Treat input as file")
    input = data.table::fread(input, header = FALSE)
  } else {
    data.table::setDT(input)
  }
  
  input = input[input[[3]] == "gene", ]
  
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  
  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  
  ensembl_id_to_genename <- data.frame(gene_id = gene_id,
                                       gene_name = gene_name,
                                       stringsAsFactors = FALSE)
  return(ensembl_id_to_genename)
}

ensembl_id_to_genename <- get_map("/Users/namuhanz/cancer/transcripts/gencode.v38lift37.annotation.gtf")

gtf_ensembl_id <- substr(ensembl_id_to_genename[,1],1,15)
ensembl_id_to_genename <- data.frame(gtf_ensembl_id,ensembl_id_to_genename[,2])
colnames(ensembl_id_to_genename) <- c("ensembl_id", "gene_id")
# write.csv(ensembl_id_to_genename, file = "../Ensembl_ID-To_Genename.csv")

#######################################################################################################
### merge data
mergeRawCounts <- merge(ensembl_id_to_genename, RawCounts, by="ensembl_id")
### order the data by gene_id
mergeRawCounts <- mergeRawCounts[order(mergeRawCounts[,"gene_id"]),]
###setup index by gene_id
index <- duplicated(mergeRawCounts$gene_id)
mergeRawCounts <- mergeRawCounts[!index,]
###take gene_name as row name
rownames(mergeRawCounts) <- mergeRawCounts[,"gene_id"]
Counts_expMatrix <- mergeRawCounts[,-c(1:2)]
###save the file
# write.csv(Counts_expMatrix, file = "../Counts_expMatrix.csv")

#######################################################################################################
###get access to the barcode
library(TCGAbiolinks)

query <- GDCquery(project = c("TCGA-LGG", "TCGA-GBM"),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts") #change every analysis

SamplesDown <- getResults(query, cols = c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = SamplesDown, typesample = "TP")
dataSmNT <- TCGAquery_SampleTypes(barcode = SamplesDown, typesample = "NT")

library(dplyr)
a <- intersect(dataSmNT, colnames(Counts_expMatrix))
b <- intersect(dataSmTP, colnames(Counts_expMatrix))

Counts <- data.frame(c(Counts_expMatrix[,a],Counts_expMatrix[,b]))
rownames(Counts) <- row.names(Counts_expMatrix)
colnames(Counts) <- c(a, b)

#######################################################################################################
###filter out and keep tumor & normal samples from the same patient
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
