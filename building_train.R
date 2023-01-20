kmers <- read.csv("output.tsv", sep = "\t", row.names = 1 )
cond <- read.csv("samples_cond", sep = ",", header = 0)

library(dplyr)
kmers <- kmers %>% select_if(~ !any(is.na(.)))

kmers <- as.data.frame(t(kmers))

kmers$class <- NA
for (i in 1:nrow(kmers)) {
  if(row.names(kmers)[i] %in% cond[1,])
    kmers[i,"class"] <- 0
  if(row.names(kmers)[i] %in% cond[2,])
    kmers[i,"class"] <- 1
}

kmers <- kmers[complete.cases(kmers), ]

write.csv(kmers, file = "train.tsv")
