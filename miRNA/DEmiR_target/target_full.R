setwd("/Users/namuhanz/cancer/miRNA/miRNA targets/file")
library(multiMiR)
library(miRNAmeConverter)

### data read: de list and sex-specific list #####
demir1 <- read.table("./c1.xls", header = T)
demir2 <- read.table("./c2.xls", header = T)
demir3 <- read.table("./c3.xls", header = T)
demir1_up <- demir1[demir1$up.down=="Up",]
demir1_down <- demir1[demir1$up.down=="Down",]
demir2_up <- demir2[demir2$up.down=="Up",]
demir2_down <- demir2[demir2$up.down=="Down",]
demir3_up <- demir3[demir3$up.down=="Up",]
demir3_down <- demir3[demir3$up.down=="Down",]

female_up_mir <- read.csv("./female_specific_up.csv", header = T)
female_down_mir <- read.csv("./female_specific_down.csv", header = T)
male_up_mir <- read.csv("./male_specific_up.csv", header = T)
male_down_mir <- read.csv("./male_specific_down.csv", header = T)

### convert miRNA name version against miRBase ############
miR_List <- list(demir1_up,demir1_down,demir2_up,demir2_down,demir3_up,demir3_down,
              female_up_mir,female_down_mir,male_up_mir,male_down_mir)
miR_target <- list(de1_up=character(),de1_down=character(),
                   de2_up=character(),de2_down=character(),
                   de3_up=character(),de3_down=character(),
                   female_up=character(),female_down=character(),
                   male_up=character(),male_down=character())

nc <- MiRNANameConverter()
for(i in 1:length(miR_List)){
  new_id <- translateMiRNAName(nc, miR_List[[i]]$gene_id, versions = c(17,20,21),sequenceFormat = 1,verbose = F)
  miR_target[[i]] <- unique(c(new_id$input,new_id$v17.0,new_id$v20.0,new_id$v21.0, miR_List[[i]]$gene_id))
}

### miRNA targets #######
targets <- list(de1up=character(), de1down=character(),
                de2up=character(), de2down=character(),
                de3up=character(), de3down=character(),
                femaleup=character(), femaledown=character(),
                maleup=character(),maledown=character())
for(j in 1:length(miR_target)){
  res <- get_multimir(org     = 'hsa',
                      mirna   = miR_target[[j]],
                      table   = 'validated',
                      summary = TRUE)
  targets[[j]] <- unique(res@data$target_symbol)
  filename = paste("target_",names(targets)[j],".csv", sep="")
  write.csv(targets[[j]], file = filename)
}

print(c(length(targets[[1]]),length(targets[[2]]),
        length(targets[[3]]),length(targets[[4]]),
        length(targets[[5]]),length(targets[[6]]),
        length(targets[[7]]),length(targets[[8]]),
        length(targets[[9]]),length(targets[[10]])
      ))

 


