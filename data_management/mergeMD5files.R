#!/usr/bin/env Rscript

args <- commandArgs(trailing=TRUE)
input_path <- args[1]
output_path <- args[2]

df1 <- read.table(input_path, sep="\t", fill=TRUE, stringsAsFactors=FALSE, header=FALSE, quote="")
colnames(df1) <- c("file", "asize", "dsize", "human_readable", "file_extension", "md5")
table1 <- table(df1$md5)
dups <- names(table1[which(table1>1)])
df1$asize <- as.numeric(df1$asize)
df2 <- df1[df1$md5 %in% dups, ]
df2 <- df2[order(df2$asize, df2$md5, decreasing = TRUE), ]
write.table(df2, output_path, sep="\t", quote=FALSE, row.names=FALSE)
