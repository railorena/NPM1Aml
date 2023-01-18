library(seqinr)

seqs <- read.csv("output.tsv", sep = ",", header = T)
nrow(seqs)
seqnames <- paste0("seq_",rownames(seqs))
seqs <- as.list(seqs$tag)
names(seqs) <- seqnames
write.fasta(sequences=seqs,names=names(seqs),file.out="output.fa")
