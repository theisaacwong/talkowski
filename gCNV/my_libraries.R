if(unname(Sys.info()["sysname"]) != "Linux"){
  library("rgl")  
  library("parallel")
  library("plyr")
  library("factoextra")
  library("clusterSim")
  library("cluster")
  library("clValid")
  # library("GenomicRanges")
  library("gridExtra")
  library("wesanderson")
  library("wesanderson")
  library("factoextra")
  library("ClusterR")
  library("dbscan")
  library("class")
 # library("ggpubr")
}

library("stringr")
library("ggplot2")
library("data.table")
library("tidyverse")
library("fortunes")
library("cowsay")


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
                           hash <<- new.env(hash = TRUE, parent = emptyenv(), size = 1000L)
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
                         populate = function(dataframe, keyCol, valCol){
                           if(is.character(keyCol) | is.character(valCol)){
                             keyCol = colnames(dataframe)[which(colnames(dataframe) == keyCol) ]
                             valCol = colnames(dataframe)[which(colnames(dataframe) == valCol) ]
                           }
                           for(i in 1:nrow(dataframe)){
                             base::assign(dataframe[i, keyCol], dataframe[i, valCol], hash)
                           }
                         },
                         populate_append = function(dataframe, keyCol, valCol){
                           if(is.character(keyCol) | is.character(valCol)){
                             keyCol = colnames(dataframe)[which(colnames(dataframe) == keyCol) ]
                             valCol = colnames(dataframe)[which(colnames(dataframe) == valCol) ]
                           }
                           for(i in 1:nrow(dataframe)){
                             key <- dataframe[i, keyCol]
                             appendant <- dataframe[i, valCol]
                             if(base::exists(key, hash)){
                               base::assign(key, base::append(unname(base::get(key, hash)) ,appendant), hash)
                             } else {
                               base::assign(key, appendant, hash)
                             }
                           }
                         },
                         head = function(){
                           keys <- utils::head(ls(hash))
                           rval <- vector("list", length(keys))
                           names(rval) <- keys
                           for(k in keys){
                             rval[[k]] <- unname(base::get(k, hash))
                           }
                           return(rval)
                         },
                         random = function(n=10){
                           keys <- sample(ls(hash), n)
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
                         clear = function() { hash <<- new.env(hash = TRUE, parent = emptyenv(), size = 1000L) }
                       ) 
)


# get
"%G%" <- function(hashMap, key) {hashMap$get(key)} 



aggregate_sample_manifests <- function(wd){
  files <- list.files(wd, pattern = "^sample_", full.names = TRUE)
  
  # Read in terra workspace manifests
  samples <- lapply(files[grepl("sample_", files)], function(x) {
    df <- read.table(x, sep="\t", header=TRUE, stringsAsFactors = FALSE, na.strings = "wefwoeifjowiejfoiwejfiow", check.names = FALSE)
    df$source_workspace <- basename(x) %>% str_replace_all("sample_|.tsv", "")
    return(df)
  }) %>% rbind.fill()
  target_columns <- c("entity", "project", "collaborator", "version", "data_type", "path", "sample", "date", "participant", "gender", "source", "id", "terra", "workspace", "md5", "name")
  target_columns_regex <- paste("(", paste(target_columns, collapse = ")|("), ")", sep = "")
  target_columns_index <- grepl(target_columns_regex, colnames(samples))
  samples0 <- samples[, target_columns_index]
  samples0[is.na(samples0)] <- "NA"
  samples0[samples0 == ""] <- "NA"
  samples0$IW_ID <- seq_down(samples0) %>% paste0("GCNV_ID_", .)
  
  return(samples0)  
}

get_count_per_gene <- function(df00, genes, phenotypes, listSampleToDiagnosis) {
  df_count <- matrix(0, nrow=length(genes)*2, ncol=(1+length(phenotypes))) %>% as.data.frame( na.strings="!@#$")
  colnames(df_count) <- c("gene", phenotypes)
  temp1 <- expand.grid(genes, c("DUP", "DEL")) 
  gene_plus_svtype <- paste0(temp1$Var1, "_", temp1$Var2)
  df_count$gene <- gene_plus_svtype
  
  map_classNameToColIndex <- HashMap$new()
  lapply(seq_along(phenotypes), function(x) {
    map_classNameToColIndex$put(phenotypes[x], (x+1))  
  })
 
  map_geneNameToIndex <- HashMap$new()
  lapply(seq_along(gene_plus_svtype), function(x) {
    map_geneNameToIndex$put(gene_plus_svtype[x], (x))  
  })
  
  df00 <- df00[df00$sample %in% names(listSampleToDiagnosis), ]
  
  for(i in seq_down(df00)){
    
    curr_phenotypes <- listSampleToDiagnosis[[df00$sample[i]]]
    
    for(curr_phenotype in curr_phenotypes){
      curr_genes <- str_split(df00$genes_strict_overlap[i], ",")[[1]]
      if(curr_genes[1] == "None") {next}
      
      curr_svtype <- df00$svtype[i]
      
      for(g in curr_genes){
        
        g_row <- map_geneNameToIndex$get(paste0(g, "_", curr_svtype))
        c_col <- map_classNameToColIndex$get(curr_phenotype)
        
        df_count[g_row, c_col] <- df_count[g_row, c_col] + 1
      }
    }
  }
  
  return(df_count)
}


getNormalizedCounts <- function(counts_path){
  counts_df <- read.table(counts_path, sep="\t", stringsAsFactors = FALSE, header=TRUE, fill=TRUE)
  
  counts_df[is.na(counts_df)] <- 0
  counts_matrix <- t(as.matrix(counts_df))
  indexes_to_remove <- rep(FALSE, nrow(counts_matrix))
  rown <- rownames(counts_matrix)
  chr <- do.call(rbind, str_split(rown, "_"))[,1]
  bool_y <- chr == "X" | chr == "chrX" | chr == "ChrX"
  bool_x <- chr == "Y" | chr == "chrY" | chr == "ChrY"
  indexes_to_remove[bool_y | bool_x] <- TRUE
  mydata_filtered <- counts_matrix[!indexes_to_remove,]
  mydataNormalized <- t(t(mydata_filtered)/colSums(mydata_filtered, na.rm = TRUE))
  
  return(mydataNormalized)
}

getCountsPCA <- function(counts_path, dimensions = 3){
  pca <- prcomp(t(getNormalizedCounts(counts_path)), rank = dimensions)
  return(pca)
}

clusterPrimaryClusters <- function(pca_loads, n_clusters=30, my_seed=123, wd ="", USE_METRICS=FALSE) {
  set.seed(my_seed)
  
  pca_loadings <- pca_loads
  n_clusters <- n_clusters
  hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid") 
  distance_methods <- c("euclidean", "maximum", "canberra")
  mink_p <- 10
  results <- expand.grid(hclust_methods, distance_methods)
  colnames(results) <- c("agglomeration", "distance")
  results$db <- 0
  results$silhouette <- 0
  results$dunn <- 0
  
  dist_eucl <- dist(pca_loadings[,1:3], method="euclidean") 
  print("eucl!")
  dist_maxi <- dist(pca_loadings[,1:3], method="maximum") 
  print("maxi!")
  dist_canb <- dist(pca_loadings[,1:3], method="canberra") 
  print("canb!")  
  # dist_mink <- dist(pca_loadings[,1:3], method="minkowski", p=mink_p)
  # print("mink!")


  cuts <- lapply(seq_down(results), function(i) {
    print(i)
    dist_mat <- switch(results$distance[i] %>% as.character(), "euclidean"=dist_eucl, "maximum"=dist_maxi, "canberra"=dist_canb)
    hclust <- hclust(dist_mat, method=results$agglomeration[i])
    cut_avg <- cutree(hclust, k=n_clusters)
    temp_names_srt <- names(srt(cut_avg))
    return(sapply(unname(cut_avg), function(x) which(x == temp_names_srt)))
  })
  
  # cuts <- vector("list", n_clusters)  
  # for(i in 1:nrow(results)){
  #   print(i)
  #   dist_mat <- switch(results$distance[i] %>% as.character(), "euclidean"=dist_eucl, "maximum"=dist_maxi, "canberra"=dist_canb)
  #   hclust <- hclust(dist_mat, method=results$agglomeration[i]) 
  #   cut_avg <- cutree(hclust, k=n_clusters)
  #   temp_names_srt <- names(srt(cut_avg))
  #   cuts[[i]] <- sapply(unname(cut_avg), function(x) which(x == temp_names_srt))
  #   
  #   if(USE_METRICS){
  #     results$db[i] <- index.DB(pca_loadings[,1:3], cut_avg, centrotypes="centroids", p=3)$DB
  #     results$silhouette[i] <- mean(silhouette(cut_avg, dist_mat)[,3])
  #     results$dunn[i] <- dunn(dist_mat, cut_avg)    
  #   }
  # }
  # 
  # if(USE_METRICS){
  #   results <- results[order(results$silhouette, decreasing = TRUE),]
  #   rval <- results
  #   rval$db <- order(results$db, decreasing = FALSE)
  #   rval$dunn <- order(results$dunn, decreasing = TRUE)
  #   rval$silhouette <- order(results$silhouette, decreasing = TRUE)
  #   rval$sum <- rval$db + rval$dunn + rval$silhouette
  #   rval <- rval[order(rval$sum),]
  #   write.table(rval , paste0(wd, "clustering_metrics", n_clusters, ".tsv"), sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  # }
  
  df_cuts <- do.call(cbind, cuts)
  ct <- as.data.frame(do.call(rbind, lapply(cuts, table)))
  ct <- ct[order(ct[,1]), ]
  r_val <- list(df_cuts, ct)
  names(r_val) <- c("cluster_labels", "cluster_counts")
  return(r_val)
  
}


relabelSmallClusters <- function(pca, n_clusters, sorted_clusters, seed = 123, threshold = 200){
  set.seed(123)
  threshold <- 200
  cluster_centers <- do.call(rbind, lapply(1:n_clusters, function(i){
    message(i)
    colMeans(pca$x[sorted_clusters==i, ])  
  }) ) %>% as.data.frame()
  
  tab <- table(sorted_clusters)
  which_clusters_to_relabel <- which(tab < threshold) %>% unname
  clusters_to_relabel <- pca$x[sorted_clusters %in% which_clusters_to_relabel, ]
  bool_cluster_relabel <- sorted_clusters %in% which_clusters_to_relabel
  main_cluster_centers <- cluster_centers[-c(which_clusters_to_relabel),]
  
  
  nearest_clusters <- lapply(1:nrow(pca$x), function(i){
    if(bool_cluster_relabel[i]){ 
      dists <- as.matrix(dist(rbind(pca$x[i, ], main_cluster_centers)))[1, -1]
      return(names(sort(dists))[1] %>% as.integer())
    } else {
      return(sorted_clusters[i])
    }
  }) %>% unlist %>% as.numeric()
  
  
  cohort_size <- 200
  memb <- do.call(rbind, lapply(sort(unique(nearest_clusters)), function(i){
    if(length(which(sorted_clusters==i)) <= cohort_size){
      df <- data.frame(`membership:sample_set_id`=paste0("cluster_super_", i, "_COHORT"), 
                       samples = row.names(pca$x)[which(nearest_clusters==i)], 
                       check.names = FALSE)
      return(df)
    } else {
      cohort_subset_indexes <- base::sample(which(sorted_clusters==i), cohort_size, replace = FALSE)
      nearest_points_in_cluster_indexes <- which(nearest_clusters==i)
      case_subset_indexes <- base::setdiff(nearest_points_in_cluster_indexes, cohort_subset_indexes)
      
      cohort_label <- paste0("cluster_super_", i, "_COHORT")
      case_label <- paste0("cluster_super_", i, "_CASE")
      
      cohort_samples <- row.names(pca$x)[cohort_subset_indexes]
      case_samples <- row.names(pca$x)[case_subset_indexes]
      
      df <- data.frame(`membership:sample_set_id`=c(
        rep(cohort_label, length(cohort_samples)), 
        rep(case_label, length(case_samples)) ), 
        samples=c(
          cohort_samples,
          case_samples),
        stringsAsFactors = FALSE, 
        check.names = FALSE)
      return(df)
    }
  }) )
  
  r_val <- list(memb, nearest_clusters)
  names(r_val) <- c("membership_file", "merged_clusters")
  return(r_val)
}


clusterOnSex <- function(wd, membership_file, seed = 123){
  set.seed(seed)
  files <- list.files(wd, pattern=".exons.counts.tsv", full.names = TRUE, recursive = TRUE)
  
  cluster <- makeCluster(detectCores()-1)
  clusterExport(cl=cluster, varlist=c("files", "str_replace_all"))
  counts_dfs <- parLapply(cluster, files, function(x){
    cbind(read.table(x, sep="\t", header=TRUE, stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "count"), comment.char = "@"),
          str_replace_all(basename(x), ".exons.counts.tsv", ""))
  })
  stopCluster(cluster)
  
  
  names(counts_dfs) <- str_replace_all(basename(files), ".exons.counts.tsv", "")
  
  total_counts <- lapply(counts_dfs, function(x){
    x <- x[!(x$start<=2781479 & x$end>=10001),] # remove PAR1
    x <- x[!((x$chr=="chrX" & x$start<=156030895 & x$end>=155701383) | 
               (x$chr=="chrY" & x$start<=57217415 & x$end>=56887903)),] # remove PAR2
    data.frame(xcount = sum(x$count[x$chr=="chrX"]), 
               ycount = sum(x$count[x$chr=="chrY"]), 
               sample = as.character(x[1,5]), stringsAsFactors = FALSE)
  })
  
  merged_counts <- do.call(rbind, total_counts)
  merged_counts$x <- merged_counts$xcount / median(merged_counts$xcount)
  merged_counts$y <- merged_counts$ycount / median(merged_counts$ycount)
  
  memb <- read.table(membership_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  memb <- memb[str_detect(memb$`membership:sample_set_id`, "cluster_super_.*[COHORT_CASE]"), ]
  
  sampleToCluster <- HashMap$new()
  sampleToCluster$populate(memb, 2, 1)
  
  merged_counts$membership <- sapply(merged_counts$sample, function(x){
    sampleToCluster$get(x)
  })
  
  merged_counts$super_cluster <- merged_counts$membership %>% str_replace_all("_(COHORT|CASE)", "") %>% str_replace_all("_super", "")
  
  g1 <- ggplot(merged_counts, aes(x=x, y=y, color=super_cluster)) + 
    geom_point() + 
    xlab("normalized chrX counts") + 
    ylab("normalized chrY counts") + 
    ggtitle("normalized chrX and chrY counts")
  print(g1)
  
  
  n_super_clusters <- merged_counts$super_cluster %>% unique %>% length
  super_clusters <- merged_counts$super_cluster %>% unique
  merged_counts$sex_cluster <- "-1"
  merged_counts$sex_cluster_sub <- "-1"
  
  clustering_method <- "hclust"
  
  for(i in 1:n_super_clusters){
    n_sub_cluster <- 2
    current_super_cluster <- merged_counts[merged_counts$super_cluster == super_clusters[i],]
    
    if(clustering_method == "hclust"){ # hclust is generally better using ward.D2 99% of the time so is ahrd coded here
      hclust <- hclust(dist(current_super_cluster[, c("x", "y")]), method="ward.D2")
      cut_avg <- cutree(hclust, k=n_sub_cluster)
      labels <- c("A", "B")[cut_avg]
      sex_labels <- labels
    } else if(clustering_method == "GMM"){
      gclust <- GMM(current_super_cluster[, c("x", "y")], n_sub_cluster, km_iter = 20, em_iter=200)
      pr <- predict_GMM(current_super_cluster[, c("x", "y")], gclust$centroids, gclust$covariance_matrices, gclust$weights)
      labels <- c("A", "B")[1 + pr$cluster_labels]
      sex_labels <- labels
    } else if(clustering_method == "DBSCAN"){
      dclust <- dbscan(current_super_cluster[, c("x", "y")], eps=0.1, minPts = 10); table(dclust$cluster); lab <- dclust$cluster + 1
      n_clusters <- length(unique(lab))
      x <- current_super_cluster[, c("x", "y")]
      sorted_clusters <- sapply(unname(lab), function(x) which(x == names(srt(lab))))
      cluster_centers <- do.call(rbind, lapply(1:n_clusters, function(i){
        if(nrow(x[sorted_clusters==i, ]) == 1){
          return(x[sorted_clusters==i, ])    
        } else {
          colMeans(x[sorted_clusters==i, ])  
        }
      }) ) %>% as.data.frame()
      
      tab <- table(sorted_clusters)
      which_clusters_to_relabel <- sorted_clusters[!sorted_clusters %in% c(1,2)] %>% unique()
      bool_cluster_relabel <- sorted_clusters %in% which_clusters_to_relabel
      main_cluster_centers <- cluster_centers[-c(which_clusters_to_relabel),]
      
      nearest_clusters <- lapply(1:nrow(x), function(i){
        if(!bool_cluster_relabel[i]){ 
          return(sorted_clusters[i])
        } else {
          dists <- as.matrix(dist(rbind(x[i, ], main_cluster_centers)))[1, -1]
          return(as.numeric(names(dists)[which(dists==min(dists))][1]))
        }
      }) %>% unlist %>% as.numeric()
      
      labels <- c("A", "B")[nearest_clusters]
      sex_labels <- labels
    }
    
    
    # this is quite bad rn, B, C denote quality
    x_mean_A <- mean(current_super_cluster$x[labels=="A"])
    x_mean_B <- mean(current_super_cluster$x[labels=="B"])
    y_mean_A <- mean(current_super_cluster$y[labels=="A"])
    y_mean_B <- mean(current_super_cluster$y[labels=="B"])
    if(x_mean_A < x_mean_B & y_mean_A > y_mean_B){
      sex_labels[sex_labels=="A"] <- "Y"
      sex_labels[sex_labels=="B"] <- "X"
    } else if(x_mean_A > x_mean_B & y_mean_A < y_mean_B) {
      sex_labels[sex_labels=="A"] <- "X"
      sex_labels[sex_labels=="B"] <- "Y"
    } else if(y_mean_A > y_mean_B){
      sex_labels[sex_labels=="A"] <- "YB"
      sex_labels[sex_labels=="B"] <- "XB"
    } else if(y_mean_A < y_mean_B) {
      sex_labels[sex_labels=="A"] <- "XB"
      sex_labels[sex_labels=="B"] <- "YB"
    } else if(x_mean_A < x_mean_B){
      sex_labels[sex_labels=="A"] <- "YC"
      sex_labels[sex_labels=="B"] <- "XC"
    } else if(x_mean_A > x_mean_B) {
      sex_labels[sex_labels=="A"] <- "XC"
      sex_labels[sex_labels=="B"] <- "YC"
    } 
    sub_labels <- sex_labels
    labels <- sub_labels
    current_super_cluster$sub_labels <- sex_labels

    g2 <- ggplot(current_super_cluster, aes(x=x, y=y, color=sub_labels)) + 
        geom_point() + 
        xlab("normalized chrX counts") + 
        ylab("normalized chrY counts") + 
        ggtitle(super_clusters[i]) 
    print(g2)
        
    
    merged_counts$sex_cluster[merged_counts$super_cluster == super_clusters[i]] <- paste0(super_clusters[i], "_", labels) 
    
    for(k in unique(labels)){
      if(length(which(labels==k)) <= 200){
        sub_labels[which(labels==k)] <- paste0(sub_labels[which(labels==k)], "_COHORT")
      } else {
        cohort_subset_indexes <- base::sample(which(labels==k), 200, replace = FALSE)
        case_subset_indexes <- base::setdiff(which(labels==k), cohort_subset_indexes)
        
        sub_labels[cohort_subset_indexes] <- paste0(sub_labels[cohort_subset_indexes], "_COHORT")
        sub_labels[case_subset_indexes] <- paste0(sub_labels[case_subset_indexes], "_CASE")
        
      }
    }
    
    merged_counts$sex_cluster_sub[merged_counts$super_cluster == super_clusters[i]] <- paste0(super_clusters[i], "_", sub_labels) 
  }
  
  table(merged_counts$sex_cluster_sub) %>% print
  memb <- data.frame(membership.sample_set_id=merged_counts$sex_cluster_sub, sample=merged_counts$sample)
  return(memb)
  #write.table(memb, "C:/Users/iwong/Documents/MGH/membership.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
}


fortune() %>% paste(collapse = "\n\t-") %>% say(by="random")


