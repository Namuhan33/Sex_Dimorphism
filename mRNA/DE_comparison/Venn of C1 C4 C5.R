setwd("/Users/namuhanz/cancer/transcripts/compare DE/file")
c1 <- read.table('c1.xls', header = T)
c4 <- read.table('c4.xls', header = T)
c5 <- read.table('DEG_FDR&FC_m.xls', header = T)
c1_gene <- c1$gene_id
c4_gene <- c4$gene_id
c5_gene <- c5$gene_id

library(VennDiagram)
venn.diagram(
  x = list(c1_gene, c4_gene, c5_gene),
  category.names = c(paste("~type", "(",length(c1_gene), ")", sep = ""), paste("~sex+type", "(",length(c4_gene),")", sep = ""), paste("~sex*type", "(",length(c5_gene),")", sep = "")),
  filename = 'Sex Effect Varification.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 520 , 
  width = 520 , 
  resolution = 400,
  compression = "lzw",
  # Circles
  lwd = 1,
  lty = 1,
  # col = c(alpha("#440154ff",10), alpha('#21908dff',10), alpha('#fde725ff',10)),
  col = c("darkorange3", '#21908dff', 'tomato'),
  fill = c(alpha("darkorange3",0.5), alpha('#21908dff',0.5), alpha('tomato',0.5)),
  # fill = myCol,
  # Numbers
  cex = 0.3,
  fontface = "plain",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.25,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.055, 0.075, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("darkorange3", '#21908dff', 'tomato'),
  rotation = 1,
  main = NULL,
  main.cex = 0.4,
  main.fontfamily	= "sans",
  main.pos = c(0.5,0.9),
  margin = 0.05
)
