#!/usr/bin/env Rscript
library(stringr)
args = commandArgs(trailingOnly=TRUE)

folderPath <- args[1]
files <- list.files(path=folderPath)
outFile <- args[2]

data <- do.call(rbind, lapply(files, function(x) cbind(x, read.table(paste(folderPath, x, sep=""), header=TRUE, comment.char = "@", stringsAsFactors = FALSE, sep = "\t", colClasses = c(rep("NULL", 3), "integer")))))
mat <- matrix(data$COUNT, ncol = length(files), byrow = TRUE)
individuals <- str_extract(files, "[^.]*")
indOne <- read.table(paste(folderPath, files[1], sep=""), header = TRUE, comment.char = "@", stringsAsFactors = FALSE, sep="\t")
coords <- paste(indOne$CONTIG,":", indOne$START, "-",indOne$END, sep="")
mat <- t(mat)

row.names(mat) <- paste("HG", individuals, sep="")
colnames(mat) <- coords

# everything is ok until the next line

write.table(mat, outFile, col.names = NA, quote = FALSE, sep="\t")