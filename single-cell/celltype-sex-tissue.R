x1 <- read.csv("../metabeforenorm.csv", header = T)
x1 <- x1[,-1]
x_tcell <- x1[celltype=='Tcell',]
x_bcell <- x1[celltype=='Bcell',]
x_epi <- x1[celltype=='Epithelial',]
x_endo <- x1[celltype=='Endothelial',]
x_mac <- x1[celltype=='Macrophage',]
x_fib <- x1[celltype=='Fibroblast',]
x_mast <- x1[celltype=='MastCell',]
x_na <- x1[is.na(x1$celltype),]

table(x_mast$sex.tissue)
