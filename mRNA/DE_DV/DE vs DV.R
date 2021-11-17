setwd("/Users/namuhanz/cancer/transcripts/DE vs DV/file")

dv <- read.csv("./DV_MAD.csv", header = T)
de <- read.table("./DE_FDR&FC.xls", header = T)

dv_n <- dv[dv$type=="N",]$DVgene
dv_t <- dv[dv$type=="T",]$DVgene
de_up <- de[de$up.down=="Up",]$gene_id
de_down <- de[de$up.down=="Down",]$gene_id

up_n <- intersect(dv_n,de_up)
up_t <- intersect(dv_t,de_up)
down_n <- intersect(dv_n,de_down)
down_t <- intersect(dv_t,de_down)

write.csv(up_n, file = "./DEup-DVnormal.csv")
write.csv(up_t, file = "./DEup-DVtumor.csv")
write.csv(down_n, file = "./DEdown-DVnormal.csv")
write.csv(down_t, file = "./DEdown-DVtumor.csv")

n <- length(up_n)+length(down_n)
t <- length(up_t)+length(down_t)
u <- length(up_n)+length(up_t)
d <- length(down_n)+length(down_t)

print(c(length(up_n), length(up_t),u))
print(c(length(down_n),length(down_t),d))
print(c(n,t))



library(Exact)

a <- as.integer(length(up_n))
b <- as.integer(length(up_t))
c <- as.integer(length(down_n))
d <- as.integer(length(down_t))

normal <- c(a,c)
tumor <- c(b,d)
df <- cbind(normal,tumor)
rownames(df) <- c("up","down")

exact.test(df, alternative = "two.sided", npNumbers = 100,  beta = 0.001,
           method = "Z-pooled", model = "Binomial")


