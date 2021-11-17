# temp1 #########
f_u <- read.csv("./femaledm_trueup.csv", header = T)
f_u <- f_u$x
f_d <- read.csv("./femaledm_truedown.csv", header = T)
f_d <- f_d$x
m_u <- read.csv("./maledm_trueup.csv", header = T)
m_u <- m_u$x
m_d <- read.csv("./maledm_truedown.csv", header = T)
m_d <- m_d$x

female <- c(f_u,f_d)
male <- c(m_u,m_d)

write.csv(female, file = "./female_bladder_dmtt.csv")
write.csv(male, file = "./male_bladder_dmtt.csv")






# temp2 #########
new_id <- translateMiRNAName(nc, 
                             "hsa-mir-1307", 
                             versions = c(17,20,21),
                             sequenceFormat = 1,
                             verbose = F)

res <- get_multimir(org     = 'hsa',
                    mirna   = unique(c("hsa-mir-1307", new_id$input,new_id$v17.0,new_id$v20.0,new_id$v21.0)),
                    table   = 'validated',
                    summary = TRUE)

targets <- unique(res@data$target_symbol)
a <- intersect("ELAVL3", targets)






