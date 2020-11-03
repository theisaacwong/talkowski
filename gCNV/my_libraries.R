if(unname(Sys.info()["sysname"]) != "Linux"){
  library("rgl")  
  library("parallel")
  library("plyr")
  library("factoextra")
  library("clusterSim")
  library("cluster")
  library("clValid")
  library("GenomicRanges")
  library("gridExtra")
  library("wesanderson")
  library("wesanderson")
  library("factoextra")
  library("ClusterR")
  library("dbscan")
 # library("ggpubr")
}

library("stringr")
library("ggplot2")
library("data.table")
library("tidyverse")


srt <- function(x){sort(table(x), decreasing = TRUE)}

read <- function(wd, file, head=TRUE){read.table(paste0(wd, file), sep="\t", header = head, stringsAsFactors = FALSE, na.strings = "!@#$%^")}

u <- function(x) {unique(x)}
lu <- function(x) {length(unique(x))}

heda <- function(x) {head(x)}
heaD <- function(x) {head(x)}

seq_down <- function(x){ 1:nrow(x)}
lo <- function(x){sapply(x, length)==1}
ld <- function(x){sapply(x, length) %>% table}

makeFiles <- function(name, cohort, mat, gbucket){
  output_name <- paste0(cohort, "-", name, ".txt")
  colind <- which(colnames(mat)==name)
  tmp <- mat[,colind]
  tmp <- str_replace_all(tmp, '\"', "")
  files <- unlist(str_split(str_replace_all(tmp, "(\\[)|]", ""), ","))
  write.table(files, quote=FALSE, sep="\t", col.names=F, row.names=F, file=paste0("~/downloads/", output_name))
  system2("gsutil", paste0("cp ~/downloads/", output_name, " ", gbucket, output_name))
  #file.remove(paste0("~/downloads/", output_name))
}

wrapper <- function(x, ...) {
  paste(strwrap(x, ...), collapse = "\n")
}

plot_ploidy <- function(df, i=1){
  ggplot(df, aes(x=START, y=LINEAR_COPY_RATIO)) + 
    geom_point(alpha=0.2) + 
    geom_hline(aes(yintercept=2)) +
    geom_hline(aes(yintercept=1)) +
    ylim(-0.5, 4) + 
    labs(x="Start Coordinate", 
         y="linear copy ratio", 
         title=sprintf("%s; %s; ploidy: %s", 
                       df$sample[1], 
                       df$CONTIG[i], 
                       ploidy_df[ploidy_df$CONTIG==df$CONTIG[i] & ploidy_df$sample==df$sample[1],]$PLOIDY[1]) %>% wrapper(30) )
}

# initialize with
# myMap <- HashMap$new()
HashMap <- setRefClass("HashMap", 
                       fields = list(hash = "environment"), 
                       methods = list(
                         initialize = function(...) {
                           hash <<- new.env(hash = TRUE, parent = emptyenv(), size = 100L)
                           callSuper(...)
                         },
                         get = function(key) { unname(base::get(key, hash)) },
                         put = function(key, value) { base::assign(key, value, hash) },
                         append = function(key, appendant){
                           if(base::exists(key, hash)){
                             base::assign(key, base::append(unname(base::get(key, hash)) ,appendant), hash)
                           } else {
                             base::assign(key, appendant, hash)
                           }
                         },
                         containsKey = function(key) { base::exists(key, hash) },
                         keySet = function() { base::ls(hash) },
                         toString = function() {
                           keys <- ls(hash)
                           rval <- vector("list", length(keys))
                           names(rval) <- keys
                           for(k in keys){
                             rval[[k]] <- unname(base::get(k, hash))
                           }
                           return(rval)
                         },
                         size = function() { length(hash) },
                         isEmpty = function() { length(hash)==0 },
                         remove = function(key) { .self$hash[[key]] <- NULL }, 
                         clear = function() { hash <<- new.env(hash = TRUE, parent = emptyenv(), size = 100L) }
                       ) 
)


# get
"%G%" <- function(hashMap, key) {hashMap$get(key)} 
