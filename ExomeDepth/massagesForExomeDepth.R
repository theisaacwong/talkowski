#!/usr/bin/env Rscript
library('stringr')
library('ExomeDepth')

args = commandArgs(trailingOnly=TRUE)
args <- c("/data/talkowski/iwong/CMG/countmat/", "/data/talkowski/iwong/CMG/ExomeDepth/ExomeDepth_output.txt")
inputFolderPath <- args[1]
outputFilePath <- args[2]

files <- list.files(path=inputFolderPath)
indOne <- read.table(paste(inputFolderPath, files[1], sep=""), header = TRUE, comment.char = "@", stringsAsFactors = FALSE, sep="\t")

space <- indOne$CONTIG
space[space=='X'] <- 23
space[space=='Y'] <- 24
start <- indOne$START
end <- indOne$END
width <- end - start + 1
names <- paste(indOne$CONTIG, indOne$START, sep="_")
GC <- rep(NA, length(space))
chromosome <- space

exomes <- matrix(nrow = length(indOne$COUNT), ncol = length(files))
for(i in 1:length(files)){
  currInd <- read.table(paste(inputFolderPath, files[i], sep=""), header = TRUE, comment.char = "@", stringsAsFactors = FALSE, sep="\t")
  exomes[,i] <- currInd$COUNT
}

Exomes <- data.frame(exomes)
colnames(Exomes) <- paste(rep("Exome", length(files)), 1:(length(files)), sep="")

df <- data.frame(space, start, end, width, names, GC, Exomes, chromosome)


test <- new('ExomeDepth', test=df$Exome2, reference=df$Exome3, formula='cbind(test, reference) ~ 1', subset.for.speed=seq(1, nrow(df), 100))
show(test)

my.test <- df$Exome1
my.ref.samples <- colnames(Exomes)[-1]
my.references.set <- as.matrix(df[, my.ref.samples])
my.choice <- select.reference.set(test.counts=my.test, reference.counts = my.references.set, bin.length = (df$end -df$start)/1000, n.bins.reduced = 10000 )
print(my.choice[[1]])

my.matrix <- as.matrix(df[,my.choice$reference.choice, drop=FALSE])
my.reference.selected <- apply(X = my.matrix, MAR = 1, FUN = sum)

all.exons <- new('ExomeDepth', test=my.test, reference=my.reference.selected, formula='cbind(test, reference) ~ 1')


all.exons <- CallCNVs(x=all.exons, transition.probability = 10^-4, chromosome = df$space, start=df$start, end=df$end, name=df$names)

head(all.exons@CNV.calls)

write.table(file = outputFilePath, x = all.exons@CNV.calls,row.names = FALSE, sep="\t")



