setwd("/Users/namuhanz/cancer/methylation/CG classification")
fulllist <- read.csv("jhu-usc.edu_KICH.HumanMethylation450.1.lvl-3.TCGA-KO-8415-01A-11D-2312-05.gdc_hg38.txt", sep = "\t", header = T)
fulllist <- fulllist[fulllist$Gene_Symbol != ".",]
cpgnotisland <- fulllist[fulllist$Feature_Type != "Island",]
cpgisland <- fulllist[fulllist$Feature_Type == "Island",]
cpgisland <- cpgisland[,c(-2,-3,-4,-5,-7,-8,-11)]
expandedcpgisland <- data.frame(matrix(NA, nrow = 1, ncol = ncol(cpgisland)+1))
colnames(expandedcpgisland) <- c(colnames(cpgisland),"position")
expandedcpgisland <- expandedcpgisland[,-3]

for(i in 1:nrow(cpgisland)){ # read in every row of cgid
  genesymbols <- unlist(strsplit(cpgisland[i,"Gene_Symbol"], ";"))
  position_to_TSS <- unlist(strsplit(cpgisland[i,"Position_to_TSS"], ";"))
  position <- character()
  for (j in 1:length(position_to_TSS)) {
    if(position_to_TSS[j] >= 0) {
      position[j] <- "E"
    } else {
      position[j] <- "P"
    }
  } # store the position of the cg in each of the transcripts 
  gene <- unique(genesymbols) # output = geneA, geneB, geneC
  newgenevector <- character()
  eee <- character()
  for(k in gene){
    ppp <- position[which(genesymbols==k)]
    if (length(unique(ppp)) == 1){
      numofnewgene <- 1
    } else {
      numofnewgene <- 2
    }
    newgene <- rep(k, numofnewgene)
    newgenevector <- append(newgenevector,newgene)
    eee <- append(eee,unique(ppp))
  }  
  
  for (l in 1:length(newgenevector)) {
    cgid <- cpgisland[i,"Composite.Element.REF"]
    symbol <- newgenevector[l]
    cgi <- cpgisland[i,"CGI_Coordinate"]
    posi <- eee[l]
    expanded <- c(cgid, symbol,cgi, posi)
    expandedcpgisland <- rbind(expandedcpgisland, expanded)
  }
  
  
}
















