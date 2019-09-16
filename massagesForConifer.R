#!/usr/bin/env Rscript

library("stringr")

args = commandArgs(trailingOnly=TRUE)

inputFolderPath <- args[1]
outFolderPath <- args[2]

inputFolderPath <- 'C:/Users/Isaac/OneDrive - WEOCD - Wong Family/MGH/REU_Isaac2/'
outFolderPath <- 'C:/Users/Isaac/OneDrive - WEOCD - Wong Family/MGH/temp2/'

inputFolderPath <- "/data/talkowski/iwong/CMG/countmat/"
outFolderPath <- "/data/talkowski/iwong/CMG/conifer/data_files_02/"

if(dir.exists(outFolderPath) == FALSE){
  dir.create(outFolderPath)
  dir.create(paste(outFolderPath, "/", "RPKM_data", sep=""))
} else {
  print("error! folder already exists")
  quit()
}

files <- list.files(path=inputFolderPath)

indOne <- read.table(paste(inputFolderPath, files[1], sep=""), header = TRUE, comment.char = "@", stringsAsFactors = FALSE, sep="\t")

probes <- data.frame(chr=indOne$CONTIG, start=indOne$START, end=indOne$END, name=paste(indOne$CONTIG, indOne$START, sep="."))
write.table(probes, paste(outFolderPath, "/probes.txt", sep=""), row.names = FALSE, quote=FALSE, sep="\t")

lengths <- indOne$END - indOne$START
names <- str_extract(files, "[^.]*")
names <- str_replace_all(str_replace_all(files, "SSC-SFARI__", ""), ".tsv", "")
for(i in 1:length(files)){
  currInd <- read.table(paste(inputFolderPath, files[i], sep=""), header = TRUE, comment.char = "@", stringsAsFactors = FALSE, sep="\t")
  reads <- currInd$COUNT
  totalReads <- sum(reads)
  rpkm <- reads / (lengths/1000 * totalReads/1000000)
  df <- data.frame(reads=reads, rpkm=rpkm)
  write.table(df, paste(outFolderPath, "/", "RPKM_data/",names[i],".txt",sep=""), quote=FALSE, col.names = FALSE, sep="\t")
}