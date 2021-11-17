setwd("/Users/namuhanz/cancer/transcripts/DE vs DV/file")
library(Exact)


a <- as.integer(10)
b <- as.integer(2854)
c <- as.integer(1944)
d <- as.integer(1)

normal <- c(a,c)
tumor <- c(b,d)
df <- cbind(normal,tumor)
rownames(df) <- c("up","down")

exact.test(df, alternative = "two.sided", npNumbers = 100,  beta = 0.001,
           method = "Z-pooled", model = "Binomial")
















