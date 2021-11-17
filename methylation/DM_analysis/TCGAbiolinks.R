### tutorial: https://www.bioconductor.org/packages/devel/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html#Epigenetic_analysis
### tutorial: https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/casestudy.html#Case_study_n_3:_Integration_of_methylation_and_expression_for_ACC
###########################################################################################
setwd("/Users/namuhanz/cancer/methylation/DM analysis/file")
samplenames <- read.csv("./samplenames.csv", header = T)
samples <- samplenames[,1]

library(TCGAbiolinks)
library(SummarizedExperiment)

query.met <- GDCquery(project = c("TCGA-THCA", "TCGA-DLBC"), 
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      barcode = samples,
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)
acc.met <- GDCprepare(query = query.met,
                      save = TRUE, 
                      save.filename = "accDNAmet.rda",
                      summarizedExperiment = TRUE)
### na.omit
acc.met <- subset(acc.met,subset = (rowSums(is.na(assay(acc.met))) == 0))
### volcano plot and DM result saved in csv
acc.met3 <- TCGAanalyze_DMC(acc.met, groupCol = "sample_type",
                           group1 = "Primary Tumor",
                           group2 = "Solid Tissue Normal",
                           p.cut = 0.05,
                           diffmean.cut = 0.25,
                           legend = "State",
                           adj.method = "BH",
                           probe.names = T,
                           plot.filename = "Tumor vs Normal MetVolcano.png")
#############################################################################################
### select out hypo- and hyper- methylated CpGs separately





