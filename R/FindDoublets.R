##########################################################################
# Predict doublets based on doublet scores
##########################################################################
CallDoublets <- function(score, rate){

  cat("Predict doublets...\n", file = stderr())
  num <- floor(length(score) * rate)
  threshold <- sort(score, decreasing = T)[num]
  pred <- score > threshold
  return(which(pred))
}

##########################################################################
# Call DoubletFinder to obtain doublet scores
##########################################################################
CallDoubletFinder <- function(counts){

  cat("Initialize Seurat objects...\n", file = stderr())
  rownames(counts) <- as.character(1:dim(counts)[1])
  colnames(counts) <- as.character(1:dim(counts)[2])
  seurat <- Seurat::CreateSeuratObject(counts)
  seurat <- Seurat::NormalizeData(seurat, verbose = F)
  seurat <- Seurat::ScaleData(seurat, verbose = F)
  seurat <- Seurat::FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
  seurat <- Seurat::RunPCA(seurat, verbose = F)

  cat("Search for best k in knn...\n", file = stderr())
  sink('NUL')
  sink(stdout(), type = "message")
  sweep.vector <- DoubletFinder::paramSweep_v3(seurat, PCs = 1:10, sct = FALSE)
  sweep.table <- DoubletFinder::summarizeSweep(sweep.vector, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.table)
  sink(NULL, type="message")
  sink()

  cat("Calculate doublet scores...\n", file = stderr())
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  pK <- as.numeric(levels(pK))[pK]
  sink('NUL')
  sink(stdout(), type = "message")
  seurat <- DoubletFinder::doubletFinder_v3(seurat, PCs = 1:10, pN = 0.25, pK = pK,
                             nExp = 0.1, reuse.pANN = FALSE, sct = FALSE)
  sink(NULL, type="message")
  sink()

  score <- seurat@meta.data[[4]]
  return(score)
}

##########################################################################
# Call cxds, bcds, or hybrid to obtain doublet scores
##########################################################################
Callscds <- function(counts, method){

  rownames(counts) <- as.character(1:dim(counts)[1])
  colnames(counts) <- as.character(1:dim(counts)[2])
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count))

  sink('NUL')
  if(method=='cxds'){
    sce <- scds::cxds(sce)
    CD <- SummarizedExperiment::colData(sce)
    score <- CD$cxds_score
  }

  if(method=='bcds'){
    sce <- scds::bcds(sce)
    CD <- SummarizedExperiment::colData(sce)
    score <- CD$bcds_score
  }

  if(method=='hybrid'){
    sce <- scds::cxds_bcds_hybrid(sce)
    CD <- SummarizedExperiment::colData(sce)
    score <- CD$hybrid_score
  }
  sink()
  return(score)
}

##########################################################################
# Call Scrublet to obtain doublet scores
##########################################################################
CallScrublet <- function(counts){

  if(!reticulate::py_module_available("scrublet")){
    cat("Install Scrublet...\n", file = stderr())
    reticulate::py_install("scrublet", pip = T)
  }

  scr <- reticulate::import('scrublet', delay_load = T)
  result <- scr$Scrublet(counts_matrix = BiocGenerics::t(counts), expected_doublet_rate = 0.167)

  reticulate::py_capture_output(
    result <- result$scrub_doublets(min_counts=2, min_cells=3,
                                    min_gene_variability_pctl=85, n_prin_comps=30L)
  )
  score <- as.vector(result[[1]])
  return(score)
}

##########################################################################
# Call scDblFinder to obtain doublet scores
##########################################################################
CallscDblFinder <- function(counts){

  rownames(counts) <- as.character(1:dim(counts)[1])
  colnames(counts) <- as.character(1:dim(counts)[2])
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

  sink('NUL')
  sink(stdout(), type = "message")
  sce <- scDblFinder::scDblFinder(sce)
  score <- sce$scDblFinder.score
  sink(NULL, type="message")
  sink()

  return(score)
}

##########################################################################
# Call DoubletDetection to obtain doublet scores
##########################################################################
CallDoubletDetection <- function(counts){

  if(!reticulate::py_module_available("doubletdetection")){
    reticulate::py_install("doubletdetection", pip = T)
  }

  doubletdetection <- reticulate::import('doubletdetection')
  clf <- doubletdetection$BoostClassifier(n_iters=5L,
                                          use_phenograph=FALSE,
                                          standard_scaling=TRUE)
  fit <- clf$fit(BiocGenerics::t(count))
  score <- as.vector(fit$doublet_score())
  return(score)
}

##########################################################################
# Call doubletCells to obtain doublet scores
##########################################################################
CalldoubletCells <- function(counts){

  options(warn=-1)
  score <- scran::doubletCells(counts)
  options(warn=0)
  return(score)
}

#' Obtain doublet scores and doublets
#'
#' Call different computational doublet-detection methods to obtain doublet scores and doublets.
#' @param counts
#' @param method
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
FindDoublets <- function(counts, method, rate = NULL){

  if(method == 'DoubletFinder'){
    cat("Execute DoubletFinder...\n", file = stderr())
    score <- CallDoubletFinder(counts = counts)
  }

  if(method%in%c('cxds', 'bcds', 'hybrid')){
    cat(paste("Execute", method, "...\n"), file = stderr())
    score <- Callscds(counts = counts, method = method)
  }

  if(method == 'Scrublet'){
    cat("Execute Scrublet...\n", file = stderr())
    score <- CallScrublet(counts = counts)
  }

  if(method == 'scDblFinder'){
    cat("Execute scDblFinder...\n", file = stderr())
    score <- CallscDblFinder(counts = counts)
  }

  if(method == 'DoubletDetection'){
    cat("Execute DoubletDetection...\n", file = stderr())
    score <- CallDoubletDetection(counts = counts)
  }

  if(method == 'doubletCells'){
    cat("Execute doubletCells...\n", file = stderr())
    score <- CalldoubletCells(counts = counts)
  }

  if(!is.null(rate)){
    index.doublet <- CallDoublets(score, rate)
    list.return <- list(score=score, doublet=index.doublet)
  }else{
    list.return <- list(score=score)
  }
  return(list.return)
}








