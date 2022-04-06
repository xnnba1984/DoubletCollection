#' Hyperparameter searching
#'
#' Perform a grid search across user-specified hyperparameters.
#' Return a combination of hyperparameters that optimize AUPRC or AUROC.
#' @param count A scRNA-seq count matrix.
#' @param method A name of doublet-detection method.
#' @param label A vector of doublet annotations.
#' @param type "AUPRC" or "AUROC".
#' @param n_neighbors Number of nearest neighbors in KNN classifier (Scrublet).
#' @param min_gene_variability_pctl The top percentile of highly variable genes (Scrublet).
#' @param n_prin_comps Number of principal components used to construct KNN classifer (Scrublet).
#' @param nfeatures Number of highly variable genes (DoubletFinder).
#' @param PCs Number of principal components used to construct KNN classifer (DoubletFinder).
#' @param nf Number of highly variable genes (scDblFinder).
#' @param includePCs Index of principal components to include in the predictors (scDblFinder).
#' @param max_depth Maximum depth of decision trees (scDblFinder).
#' @param k Number of nearest neighbors in KNN classifier (doubletCells).
#' @param d Number of principal components used to construct KNN classifer (doubletCells).
#' @param ntop.cxds Number of top variance genes to consider (cxds).
#' @param ntop.bcds Number of top variance genes to consider (bcds).
#' @param n_components Number of principal components used for clustering (DoubletDetection).
#' @param n_top_var_genes Number of highest variance genes to use (DoubletDetection).
#' @param n_iters Number of fit operations from which to collect p-values (DoubletDetection).
#'
#' @return A list of two elements. The first is a vector of hyperparameter combination that maximizes "AUPRC" or "AUROC". The second is a matrix of all hyperparameter combinations and correponding "AUPRC" or "AUROC". 
#' @export
#'
#' @examples
#' data.list <- ReadData(path = ".../real_datasets")
#' count.list <- data.list$count
#' label.list <- lapply(data.list$label, FUN = function(label){
#'  ifelse(label == 'doublet', 1, 0)
#' })
#' count <- count.list$`pbmc-1A-dm`
#' label <- label.list$`pbmc-1A-dm`
#' result.parameter.scDblFinder <- FindParameters(count, label, method = 'scDblFinder', type = 'AUPRC',
#' nf=c(500, 1000, 1500, 2000, 2500),
#' includePCs=c(3, 4, 5, 6, 7),
#' max_depth=c(3, 4, 5, 6, 7))
#'
FindParameters <- function(count, label, method, type,
                           # Scrublet
                           n_neighbors=NULL, min_gene_variability_pctl=NULL, n_prin_comps=NULL,
                           # DoubletFinder
                           nfeatures=NULL, PCs=NULL,
                           # scDblFinder
                           nf=NULL, includePCs=NULL, max_depth=NULL,
                           # doubletCells
                           k=NULL, d=NULL,
                           # cxds and bcds
                           ntop.cxds=NULL, ntop.bcds=NULL,
                           # DoubletDetection
                           n_components=NULL, n_top_var_genes=NULL, n_iters=NULL
                           ){

  cat('\nSearch parameters for', method, '...\n', file = stderr())

  if(method=='Scrublet'){
    table.result <- expand.grid(n_neighbors=n_neighbors, min_gene_variability_pctl=min_gene_variability_pctl,
                          n_prin_comps=n_prin_comps)
    table.result[[type]] <- 0
    for(i in 1:nrow(table.result)){
      cat('Calculate parameter combination', i, '...\n', file = stderr())
      neighbor <- as.integer(table.result[i, 1])
      gene <- as.integer(table.result[i, 2])
      pca <- as.integer(table.result[i, 3])
      score <- CallScrublet(count, neighbor, gene, pca)
      auc <- FindAUC.Single(score, label, type); auc
      table.result[i, type] <- auc
    }
  }

  if(method=='DoubletFinder'){
    table.result <- expand.grid(nfeatures=nfeatures, PCs=PCs)
    table.result[[type]] <- 0
    for(i in 1:nrow(table.result)){
      cat('Calculate parameter combination', i, '...\n', file = stderr())
      nfeatures <- table.result[i, 1]
      PCs <- table.result[i, 2]
      score <- CallDoubletFinder(count, nfeatures, PCs)
      auc <- FindAUC.Single(score, label, type)
      table.result[i, type] <- auc
    }
  }

  if(method=='scDblFinder'){
    table.result <- expand.grid(nf=nf, includePCs=includePCs, max_depth=max_depth)
    table.result[[type]] <- 0
    for(i in 1:nrow(table.result)){
      cat('Calculate parameter combination', i, '...\n', file = stderr())
      nf <- table.result[i, 1]
      includePCs <- table.result[i, 2]
      max_depth <- table.result[i, 3]
      score <- CallscDblFinder(count, nf, includePCs, max_depth)
      auc <- FindAUC.Single(score, label, type)
      table.result[i, type] <- auc
    }
  }

  if(method=='doubletCells'){
    table.result <- expand.grid(k=k, d=d)
    table.result[[type]] <- 0
    for(i in 1:nrow(table.result)){
      cat('Calculate parameter combination', i, '...\n', file = stderr())
      f <- table.result[i, 1]
      d <- table.result[i, 2]
      score <- CalldoubletCells(count, f, d)
      auc <- FindAUC.Single(score, label, type)
      table.result[i, type] <- auc
    }
  }

  if(method=='cxds'){
    table.result <- expand.grid(ntop.cxds=ntop.cxds)
    table.result[[type]] <- 0
    for(i in 1:nrow(table.result)){
      cat('Calculate parameter combination', i, '...\n', file = stderr())
      ntop.cxds <- table.result[i, 1]
      score <- Callscds(count, method, ntop.cxds = ntop.cxds)
      auc <- FindAUC.Single(score, label, type)
      table.result[i, type] <- auc
    }
  }

  if(method=='bcds'){
    table.result <- expand.grid(ntop.bcds=ntop.bcds)
    table.result[[type]] <- 0
    for(i in 1:nrow(table.result)){
      cat('Calculate parameter combination', i, '...\n', file = stderr())
      ntop.bcds <- table.result[i, 1]
      score <- Callscds(count, method, ntop.bcds = ntop.bcds)
      auc <- FindAUC.Single(score, label, type)
      table.result[i, type] <- auc
    }
  }

  if(method=='DoubletDetection'){
    table.result <- expand.grid(n_components=n_components, n_top_var_genes=n_top_var_genes, n_iters=n_iters)
    table.result[[type]] <- 0
    for(i in 1:nrow(table.result)){
      cat('Calculate parameter combination', i, '...\n', file = stderr())
      n_components <- table.result[i, 1]
      n_top_var_genes <- table.result[i, 2]
      n_iters <- table.result[i, 3]
      score <- CallDoubletDetection(count, n_components, n_top_var_genes, n_iters)
      auc <- FindAUC.Single(score, label, type)
      table.result[i, type] <- auc
    }
  }

  index.max <- which.max(table.result[[type]])
  #return(table.result[index.max,])
  return(list(table.result[index.max,], table.result))
}


