# Speedup UCell according to https://github.com/bigomics/plaid

calculate_Uscore_fast <- function(
    matrix, features,  maxRank=1500, w_neg=1, ties.method="average",
    ncores=1, storeRanks=FALSE, force.gc=FALSE, name="_UCell"){
  
  #Make sure we have a sparse matrix
  if (!methods::is(matrix, "dgCMatrix")) {
    matrix <- Matrix::Matrix(as.matrix(matrix),sparse = TRUE)
  }
  
  #Check if all genes in signatures are present in the data matrix
  matrix <- UCell:::check_genes(matrix, features)
  
  #Do not evaluate more genes than there are
  if (!is.numeric(maxRank)) {
    stop("Rank cutoff (maxRank) must be a number")
  }
  if (maxRank > nrow(matrix)) {
    maxRank <- nrow(matrix)
  }
  
  #Weight on neg signatures must be >=0
  if (is.null(w_neg)) {w_neg <- 1}
  if (w_neg<0) {stop("Weight on negative signatures (w_neg) must be >=0")}
  
  #Signatures cannot be larger than maxRank parameter
  sign.lgt <- lapply(features, length)
  if (any(sign.lgt > maxRank)) {
    stop("One or more signatures contain more genes than maxRank parameter.
            Increase maxRank parameter or make shorter signatures")
  }
  
  
  if(ncores > 1){
    #Split into manageable chunks index
    chunk_idy <- chunk_index(seq_len(ncol(matrix)), n=ncores)
  
    if (Sys.info()['sysname'] == "Windows" | identical(.Platform$GUI, "RStudio")) {
      BPPARAM <- BiocParallel::SnowParam(workers=ncores)
    } else {
      BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)
    }
    
    meta.list <- BiocParallel::bplapply(
      X = chunk_idy, 
      BPPARAM = BPPARAM,
      FUN = function(x) {
        cells_rankings <- data_to_ranks(matrix[,x], ties.method=ties.method)
        cells_U <- u_stat_signature_list_fast(features, cells_rankings, 
                                              maxRank=maxRank, sparse=FALSE,
                                              w_neg=w_neg)
        colnames(cells_U) <- paste0(colnames(cells_U),name)
  
        if (storeRanks==TRUE){
          gene.names <- as.character(as.matrix(cells_rankings[,1]))
          #make sparse (rank=0 means rank>=maxRank)
          cells_rankings[cells_rankings>=maxRank] <- 0
          ranks.sparse <- Matrix::Matrix(as.matrix(
            cells_rankings),sparse = TRUE)
          dimnames(ranks.sparse)[[1]] <- gene.names
          if (force.gc) {
            cells_rankings <- NULL
            gc()
          }
          return(list(cells_rankings=ranks.sparse, cells_U=cells_U))
        } else {
          if (force.gc) {
            cells_rankings <- NULL
            gc()
          }
          return(list(cells_U=cells_U))
        }
        
      })
  }else{
    cells_rankings <- data_to_ranks(matrix, ties.method=ties.method)
    cells_U <- u_stat_signature_list_fast(features, cells_rankings, 
                                          maxRank=maxRank, sparse=FALSE,
                                          w_neg=w_neg)
    colnames(cells_U) <- paste0(colnames(cells_U),name)
    
    if (storeRanks==TRUE){
      gene.names <- as.character(rownames(cells_rankings))
      #make sparse (rank=0 means rank>=maxRank)
      cells_rankings[cells_rankings>=maxRank] <- 0
      ranks.sparse <- Matrix::Matrix(as.matrix(
        cells_rankings),sparse = TRUE)
      dimnames(ranks.sparse)[[1]] <- gene.names
      if (force.gc) {
        cells_rankings <- NULL
        gc()
      }
      return(list(cells_rankings=ranks.sparse, cells_U=cells_U))
    } else {
      if (force.gc) {
        cells_rankings <- NULL
        gc()
      }
      return(list(cells_U=cells_U))
    }
  }
}


u_stat_signature_list_fast <- function(sig_list, ranks_matrix, maxRank=1000,
                                       sparse=FALSE, w_neg=1) {
  if(any(grepl('[+-]$', unique(unlist(sig_list, use.names = FALSE)), perl=TRUE), na.rm = TRUE)){
    all_gene_list <- lapply(sig_list, function(x) unique(gsub('[+-]$','',x,perl=TRUE)))
    matG <- plaid::gmt2mat(all_gene_list)
    # matG <- matG[rownames(ranks_matrix),]
    
    for(i in 1:ncol(matG)) {
      sig_neg <- grep('-$', unlist(sig_list[[i]]), perl=TRUE, value=TRUE)
      sig_neg <- gsub('-$','',sig_neg,perl=TRUE)
      if (length(sig_neg)>0) 
        matG[sig_neg, i] <- -1
    }
    # gene_list <- lapply(sig_list, function(x) unique(gsub('[+-]$','',x,perl=TRUE)))
    
    u_p <- u_stat_fast(ranks_matrix, matG, maxRank=maxRank, sign=1)
    u_n <- u_stat_fast(ranks_matrix, matG, maxRank=maxRank, sign=-1)
    
    diff <- u_p - w_neg*u_n
    diff[diff<0] <- 0
    return(diff)
  }else{
    all_gene_list <- lapply(sig_list, unique)
    matG <- plaid::gmt2mat(all_gene_list)
    # matG <- matG[rownames(ranks_matrix),]
    
    u_p <- u_stat_fast(ranks_matrix, matG, maxRank=maxRank, sign=1)
    diff <- u_p
    diff[diff<0] <- 0
    return(diff)
  }
}


u_stat_fast <- function(ranks_matrix, matG, maxRank=1000, sign=1){
  # ranks_matrix <- as.matrix(ranks_matrix[,-1])
  com_genes <- intersect(rownames(matG), rownames(ranks_matrix))
  ranks_matrix <- ranks_matrix[com_genes,]
  matG <- matG[com_genes,]
  ranks_matrix[ranks_matrix>=maxRank] <- maxRank
  
  len_sig <- Matrix::colSums(matG==sign)
  rank_sum_min <- len_sig*(len_sig+1)/2
  
  matG[matG!=sign] <- 0
  matG[matG==sign] <- 1
  
  # uscore <- 1 - (Matrix:::crossprod(matG, ranks_matrix) - rank_sum_min) / (len_sig*maxRank - rank_sum_min)
  uscore <- 1 - (plaid:::chunked_crossprod(matG, ranks_matrix) - rank_sum_min) / (len_sig*maxRank - rank_sum_min)
  uscore <- t(uscore)
}


data_to_ranks <- function(matrix, ties.method="average"){
  cells_rankings <- matrixStats::colRanks(as.matrix(matrix), ties.method = ties.method, preserveShape = TRUE)
  cells_max <- matrixStats::colMaxs(cells_rankings, na.rm = TRUE) + 1
  if(length(unique(cells_max))==1){
    cells_rankings <- cells_max[1] - cells_rankings
  }else{
    cells_rankings <- t(cells_max - t(cells_rankings))
  }
  return(cells_rankings)
}


chunk_index <- function(x, n=1, start_1=TRUE, return_list=TRUE) {
  if(n>1){
    start <- ifelse(start_1, 1, min(x))
    brks <- floor(seq(1, max(x), length.out = n + 1) )
    brks[c(1, length(brks))] <- c(start-1, max(x))
    
    left  = brks[-length(brks)]+1
    right = brks[-1]
    
    if(return_list){
      # grp <- lapply(seq_len(length(left)), function(i){c(left[i], right[i])})
      grp <- lapply(seq_len(length(left)), function(i){seq(left[i], right[i], by=1)})
    }else{
      grp <- data.frame(
        left  = left,
        right = right
      )
    }
  }else{
    left  = min(x)
    right = max(x)
    
    if(return_list){
      # grp <- list(c(left, right))
      grp <- list(seq(left, right, by=1))
    }else{
      grp <- data.frame(
        left  = left,
        right = right
      )
    }
  }
  grp
}

rankings2Uscore_fast <- function(ranks_matrix, features, w_neg=1, name="_UCell") {
  #Check if all genes in signatures are present in the stored signatures
  ranks_matrix <- check_genes(ranks_matrix, features)
  
  #Weight on neg signatures must be >=0
  if (is.null(w_neg)) {w_neg <- 1}
  if (!is.numeric(w_neg) | w_neg<0) {
    stop("Weight on negative signatures (w_neg) must be >=0")}
  
  maxRank <- max(ranks_matrix)+1
  
  cells_U <- u_stat_signature_list_fast(features, ranks_matrix, 
                                        maxRank=maxRank, sparse=FALSE,
                                        w_neg=w_neg)
  colnames(cells_U) <- paste0(colnames(cells_U),name)
  
  return(list(cells_U=cells_U))
}
