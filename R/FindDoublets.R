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
CallDoubletFinder <- function(count, nfeatures=2000, PCs=10){

  rownames(count) <- as.character(1:dim(count)[1])
  colnames(count) <- as.character(1:dim(count)[2])
  seurat <- Seurat::CreateSeuratObject(count)
  seurat <- Seurat::NormalizeData(seurat, verbose = F)
  seurat <- Seurat::ScaleData(seurat, verbose = F)
  seurat <- Seurat::FindVariableFeatures(seurat, selection.method = "vst", nfeatures = nfeatures, verbose = F)
  seurat <- Seurat::RunPCA(seurat, verbose = F)

  sink('NUL')
  sink(stdout(), type = "message")

  tryCatch({
    sweep.vector <- DoubletFinder::paramSweep_v3(seurat, PCs = 1:PCs, sct = FALSE)
    sweep.table <- DoubletFinder::summarizeSweep(sweep.vector, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep.table)

    pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
    pK <- as.numeric(levels(pK))[pK]
    seurat <- DoubletFinder::doubletFinder_v3(seurat, PCs = 1:PCs, pN = 0.25, pK = pK,
                               nExp = 0.1, reuse.pANN = FALSE, sct = FALSE)
    },
    interrupt = function(e){
      sink(NULL, type="message")
      sink()
    }
  )
  sink(NULL, type="message")
  sink()
  score <- seurat@meta.data[[4]]
  return(score)
}

##########################################################################
# Call cxds, bcds, or hybrid to obtain doublet scores
##########################################################################
Callscds <- function(count, method, ntop.cxds=500, ntop.bcds=500){

  rownames(count) <- as.character(1:dim(count)[1])
  colnames(count) <- as.character(1:dim(count)[2])
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count))

  sink('NUL')

  tryCatch({
    if(method=='cxds'){
      sce <- scds::cxds(sce, ntop = ntop.cxds)
      CD <- SummarizedExperiment::colData(sce)
      score <- CD$cxds_score
      score <- as.numeric(score)  # remove vector names
    }

    if(method=='bcds'){
      sce <- scds::bcds(sce, ntop = ntop.bcds)
      CD <- SummarizedExperiment::colData(sce)
      score <- CD$bcds_score
      score <- as.numeric(score)  # remove vector names
    }

    if(method=='hybrid'){
      sce <- scds::cxds_bcds_hybrid(sce)
      CD <- SummarizedExperiment::colData(sce)
      score <- CD$hybrid_score
      score <- as.numeric(score)  # remove vector names
    }
  },
    interrupt = function(e){
      sink()
    }
  )
  sink()
  return(score)
}

##########################################################################
# Call Scrublet to obtain doublet scores
##########################################################################
CallScrublet <- function(count, n_neighbors=round(0.5*sqrt(dim(count)[2])),
                         min_gene_variability_pctl=85L, n_prin_comps=30L){

  if(!reticulate::py_module_available("scrublet")){
    cat("Install Scrublet...\n", file = stderr())
    reticulate::py_install("scrublet",pip = T)
  }

  scr <- reticulate::import('scrublet', delay_load = T)
  result <- scr$Scrublet(counts_matrix = BiocGenerics::t(count), n_neighbors=as.integer(n_neighbors))

  reticulate::py_capture_output(
    result <- result$scrub_doublets(min_counts=2, min_cells=3,
                                    min_gene_variability_pctl=as.integer(min_gene_variability_pctl),
                                    n_prin_comps=as.integer(n_prin_comps))
  )
  score <- as.vector(result[[1]])
  return(score)
}

##########################################################################
# Call scDblFinder to obtain doublet scores
##########################################################################
CallscDblFinder <- function(count, nf=1000, includePCs=5, max_depth=5){

  rownames(count) <- as.character(1:dim(count)[1])
  colnames(count) <- as.character(1:dim(count)[2])
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = count))

  sink('NUL')
  sink(stdout(), type = "message")

  tryCatch({
    sce <- scDblFinder::scDblFinder(sce, nfeatures=nf, includePCs=includePCs, max_depth=max_depth)
    score <- sce$scDblFinder.score
  },
    interrupt = function(e){
      sink(NULL, type="message")
      sink()
      print("stop")
    }
  )
  sink(NULL, type="message")
  sink()

  return(score)
}

##########################################################################
# Call DoubletDetection to obtain doublet scores
##########################################################################
CallDoubletDetection <- function(count, n_components=30, n_top_var_genes=10000, n_iters=5){
  if(!reticulate::py_module_available("doubletdetection")){
    reticulate::py_install("doubletdetection", pip = T)
  }

  doubletdetection <- reticulate::import('doubletdetection')
  clf <- doubletdetection$BoostClassifier(n_components=as.integer(n_components),
                                          n_top_var_genes=as.integer(n_top_var_genes),
                                          n_iters=as.integer(n_iters),
                                          use_phenograph=FALSE,
                                          standard_scaling=TRUE)
  fit <- clf$fit(BiocGenerics::t(count))
  score <- as.vector(fit$doublet_score())
  cat('\n')
  return(score)
}

##########################################################################
# Call doubletCells to obtain doublet scores
##########################################################################
CalldoubletCells <- function(count, k=50, d=50){
  options(warn=-1)
  score <- scran::doubletCells(count, k=k, d=d)
  options(warn=0)
  return(score)
}

#' Calculate doublet scores on single dataset
#'
#' Call different computational doublet-detection methods to calculate doublet scores on single dataset.
#' @param count A scRNA-seq count matrix.
#' @param method A name vector of doublet-detection methods.
#' @param n_neighbors The number of nearest neighbors in KNN classifier (Scrublet).
#' @param min_gene_variability_pctl The top percentile of highly variable genes (Scrublet).
#' @param n_prin_comps Number of principal components used to construct KNN classifer (Scrublet).
#' @param nfeatures Number of highly variable genes (DoubletFinder).
#' @param PCs Number of principal components used to construct KNN classifer (DoubletFinder).
#' @param nf Number of highly variable genes (scDblFinder).
#' @param includePCs The index of principal components to include in the predictors (scDblFinder).
#' @param max_depth Maximum depth of decision trees (scDblFinder).
#' @param k The number of nearest neighbors in KNN classifier (doubletCells).
#' @param d Number of principal components used to construct KNN classifer (doubletCells).
#' @param ntop.cxds Number of top variance genes to consider (cxds).
#' @param ntop.bcds Number of top variance genes to consider (bcds).
#' @param n_components Number of principal components used for clustering (DoubletDetection).
#' @param n_top_var_genes Number of highest variance genes to use (DoubletDetection).
#' @param n_iters Number of fit operations from which to collect p-values (DoubletDetection).
#'
#' @return A list of doublet scores calculated by each doublet-detection method.
#' @export
#'
#' @examples
#' score.list <- FindScores(count = count.list$`J293t-dm`, methods = c('cxds','bcds','hybrid'))
#'
FindScores <- function(count, methods,
                       # Scrublet
                       n_neighbors=round(0.5*sqrt(dim(count)[2])), min_gene_variability_pctl=85L, n_prin_comps=30L,
                       # DoubletFinder
                       nfeatures=2000, PCs=10,
                       # scDblFinder
                       nf=1000, includePCs=5, max_depth=5,
                       # doubletCells
                       k=50, d=50,
                       # cxds and bcds
                       ntop.cxds=500, ntop.bcds=500,
                       # DoubletDetection
                       n_components=30, n_top_var_genes=10000, n_iters=5
                       ){
  score.list <- list()
  for(method in methods){
    tryCatch({
      if(method == 'DoubletFinder'){
        cat("Execute DoubletFinder...\n", file = stderr())
        score <- CallDoubletFinder(count = count, nfeatures, PCs)
        score.list[[method]] <- score
      }

      if(method%in%c('cxds', 'bcds', 'hybrid')){
        cat(paste("Execute", method, "...\n"), file = stderr())
        score <- Callscds(count = count, method = method, ntop.cxds, ntop.bcds)
        score.list[[method]] <- score
      }

      if(method == 'Scrublet'){
        cat("Execute Scrublet...\n", file = stderr())
        score <- CallScrublet(count = count, n_neighbors, min_gene_variability_pctl, n_prin_comps)
        score.list[[method]] <- score
      }

      if(method == 'scDblFinder'){
        cat("Execute scDblFinder...\n", file = stderr())
        score <- CallscDblFinder(count = count, nf, includePCs, max_depth)
        score.list[[method]] <- score
      }

      if(method == 'DoubletDetection'){
        cat("Execute DoubletDetection...\n", file = stderr())
        score <- CallDoubletDetection(count = count, n_components, n_top_var_genes, n_iters)
        score.list[[method]] <- score
      }

      if(method == 'doubletCells'){
        cat("Execute doubletCells...\n", file = stderr())
        score <- CalldoubletCells(count = count, k=k, d=d)
        score.list[[method]] <- score
      }
    }, error = function(e){
      cat('\n', method, " ERROR :", conditionMessage(e), "\n")
      score.list[[method]] <<- NA
    }
    )
  }
  return(score.list)
}


#' Calculate doublet scores on multiple dataset
#'
#' Call different computational doublet-detection methods to calculate doublet scores on multiple datasets.
#' @param count.list A list of scRNA-seq count matrix.
#' @param method A name vector of doublet-detection methods.
#' @param n_neighbors The number of nearest neighbors in KNN classifier (Scrublet).
#' @param min_gene_variability_pctl The top percentile of highly variable genes (Scrublet).
#' @param n_prin_comps Number of principal components used to construct KNN classifer (Scrublet).
#' @param nfeatures Number of highly variable genes (DoubletFinder).
#' @param PCs Number of principal components used to construct KNN classifer (DoubletFinder).
#' @param nf Number of highly variable genes (scDblFinder).
#' @param includePCs The index of principal components to include in the predictors (scDblFinder).
#' @param max_depth Maximum depth of decision trees (scDblFinder).
#' @param k The number of nearest neighbors in KNN classifier (doubletCells).
#' @param d Number of principal components used to construct KNN classifer (doubletCells).
#' @param ntop.cxds Number of top variance genes to consider (cxds).
#' @param ntop.bcds Number of top variance genes to consider (bcds).
#' @param n_components Number of principal components used for clustering (DoubletDetection).
#' @param n_top_var_genes Number of highest variance genes to use (DoubletDetection).
#' @param n_iters Number of fit operations from which to collect p-values (DoubletDetection).
#'
#' @return A list of doublet scores calculated by each doublet-detection method on multiple datasets.
#' @export
#'
#' @examples
#' data.list <- ReadData(path = ".../real_datasets")
#' count.list <- data.list$count
#' methods <- c('Scrublet','doubletCells','cxds','bcds','hybrid','scDblFinder','DoubletDetection','DoubletFinder')
#' score.list.all <- FindScores.All(count.list, methods = methods)
#'
FindScores.All <- function(count.list, methods,
                           # Scrublet
                           n_neighbors=round(0.5*sqrt(dim(count)[2])), min_gene_variability_pctl=85L, n_prin_comps=30L,
                           # DoubletFinder
                           nfeatures=2000, PCs=10,
                           # scDblFinder
                           nf=1000, includePCs=5, max_depth=5,
                           # doubletCells
                           k=50, d=50,
                           # cxds and bcds
                           ntop.cxds=500, ntop.bcds=500,
                           # DoubletDetection
                           n_components=30, n_top_var_genes=10000, n_iters=5){
  score.list.all <- list()
  for(i in 1:length(count.list)){
    count <- count.list[[i]]; dim(count)
    data.name <- names(count.list)[i]; data.name
    cat('Calculating doublet scores on dataset: ', data.name, '......\n')
    score.list <- FindScores(count = count, methods = methods,
                             n_neighbors, min_gene_variability_pctl, n_prin_comps, nfeatures, PCs, nf, includePCs, max_depth,
                             k, d, ntop.cxds, ntop.bcds, n_components, n_top_var_genes, n_iters)
    score.list.all[[data.name]] <- score.list
  }
  return(score.list.all)
}


#' Call doublets on one dataset
#'
#' Call doublets based on doublet scores and a user-specified doublet rate.
#' @param score.list A list of doublet scores on one dataset.
#' @param rate A user-specified doublet rate (from 0 to 1).
#'
#' @return A list of doublet indices for different doublet-detection methods.
#' @export
#'
#' @examples
#' doublet.list <- FindDoublets(score.list = score.list, rate = .1)
#'
FindDoublets <- function(score.list, rate){

  options(warn=-1)
  doublet.list <- list()
  for(i in 1:length(score.list)){
    score <- score.list[[i]]
    method <- names(score.list)[i]
    if(is.na(score)){
      doublet.list[[method]] <- NA
    }else{
      index.doublet <- CallDoublets(score, rate)
      doublet.list[[method]] <- index.doublet
    }
  }
  options(warn=0)
  return(doublet.list)
}


#' Call doublets on multiple datasets
#'
#' Call doublets based on doublet scores and a user-specified doublet rate on multiple datasets.
#' @param score.list.all A list of doublet scores on multiple dataset.
#' @param rate A user-specified doublet rate (from 0 to 1).
#'
#' @return A list of doublet indices for different doublet-detection methods on multiple datasets.
#' @export
#'
#' @examples
#' doublet.list.all <- FindDoublets.All(score.list.all = score.list.all, rate = .1)
#'
FindDoublets.All <- function(score.list.all, rate){
  doublet.list.all <- list()
  for(i in 1:length(score.list.all)){
    score.list <- score.list.all[[i]]
    data.name <- names(score.list.all)[i]; data.name
    cat('Call doublets on dataset:', data.name, '......\n')
    doublet.list <- FindDoublets(score.list = score.list, rate = rate)
    doublet.list.all[[data.name]] <- doublet.list
  }
  return(doublet.list.all)
}


#' Call doublets on multiple datasets and doublet rates
#'
#' Call doublets based on doublet scores and user-specified doublet rates on multiple datasets.
#' @param score.list.all A list of doublet scores on multiple dataset.
#' @param rates A user-specified vector of doublet rates (from 0 to 1).
#'
#' @return A list of doublet indices for different doublet-detection methods on multiple datasets and doublet rates.
#' @export
#'
#' @examples
#' doublet.list.all.rate <- FindDoublets.All.Rate(score.list.all = score.list.all, rates = seq(0.01, 0.25, 0.01))
#'
FindDoublets.All.Rate <- function(score.list.all, rates){
  doublet.list.all.rates <- list()
  for(rate in rates){
    print(rate)
    doublet.list.all <- FindDoublets.All(score.list.all, rate)
    doublet.list.all.rates[[as.character(rate)]] <- doublet.list.all
  }
  return(doublet.list.all.rates)
}


#' Calculate AUPRC or AUROC
#'
#' Calculate AUPRC and AUROC based on doublet score and annotation on single dataset.
#' @param score A vector of doublet scores.
#' @param label A vector of doublet annotations (0/1).
#' @param type A string of "AUPRC" or "AUROC".
#'
#' @return A number of AUPRC or AUROC.
#' @export
#'
FindAUC.Single <- function(score, label, type){

  options(warn=-1)
  fg <- score[label==1]
  bg <- score[label==0]
  if(type=='AUPRC'){
    pr <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    auc <- pr$auc.integral
  }
  if(type=='AUROC'){
    roc <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    auc <- roc$auc
  }
  options(warn=0)
  return(auc)
}



#' Calculate AUPRC or AUROC for one dataset
#'
#' Calculate AUPRC and AUROC based on doublet scores and annotations for different doublet-detection methods on one dataset.
#' @param score.list A list of doublet scores on one dataset.
#' @param label A vector of 0/1 doublet annotations.
#' @param type A character of "AUPRC" or "AUROC".
#'
#' @return A list of AUPRCs or AUROCs of different doublet-detection methods.
#' @export
#'
#' @examples
#' auprc.list <- FindAUC(score.list = score.list, label = label.list$`J293t-dm`, type = 'AUPRC')
#' auroc.list <- FindAUC(score.list = score.list, label = label.list$`J293t-dm`, type = 'AUROC')
#'
FindAUC <- function(score.list, label, type){

  options(warn=-1)
  auc.list <- list()
  for(i in 1:length(score.list)){
    score <- score.list[[i]]
    method <- names(score.list)[i]
    if(is.na(score)){
      auc.list[[method]] <- NA
    }else{
      fg <- score[label==1]
      bg <- score[label==0]
      if(type=='AUPRC'){
        pr <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
        auprc <- pr$auc.integral
        auc.list[[method]] <- auprc
      }
      if(type=='AUROC'){
        roc <- PRROC::roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
        auroc <- roc$auc
        auc.list[[method]] <- auroc
      }
    }
  }
  options(warn=0)
  return(auc.list)
}


#' Calculate AUPRC or AUROC for mutiple datasets
#'
#' Calculate AUPRC and AUROC based on doublet scores and annotations for different doublet-detection methods on multiple datasets.
#' @param score.list.all A list of doublet scores on multiple datasets.
#' @param label.list A list of vectors of 0/1 doublet annotations.
#' @param type A character of "AUPRC" or "AUROC".
#'
#' @return A list of AUPRCs or AUROCs of different doublet-detection methods on multiple datasets.
#' @export
#'
#' @examples
#' auprc.list.all <- FindAUC.All(score.list.all, label.list, 'AUPRC')
#'
FindAUC.All <- function(score.list.all, label.list, type){
  auc.list.all <- list()
  for(i in 1:length(score.list.all)){
    score.list <- score.list.all[[i]]
    label <- label.list[[i]]; table(label)
    data.name <- names(score.list.all)[i]; data.name
    cat('Calculating', type, 'on dataset:', data.name, '......\n')
    auc.list <- FindAUC(score.list = score.list, label = label, type = type)
    auc.list.all[[data.name]] <- auc.list
  }
  return(auc.list.all)
}


#' Calculate precision, recall, or true negative rate (TNR) of identified doublets
#'
#' Calculate precision, recall, or TNR based on identified doublet indices and true doublet annotations.
#' @param doublet.list A list of doublet indices for different doublet-detection methods.
#' @param label A vector of 0/1 doublet annotations.
#' @param type A character of "precision", "recall", or "TNR".
#'
#' @return A list of precision, recall, or TNR for different doublet-detection methods.
#' @export
#'
#' @examples
#' acc.list <- FindACC(doublet.list = doublet.list, label = label.list$`J293t-dm`, type = 'precision')
#'
FindACC <- function(doublet.list, label, type){

  acc.list <- list()

  for(i in 1:length(doublet.list)){

    size <- length(label)
    pred.doublet <- doublet.list[[i]]
    pred.singlet <- setdiff(1:size, pred.doublet)
    truth.doublet <- which(label==1)
    truth.singlet <- which(label==0)
    method <- names(doublet.list)[i]
    tp <- intersect(pred.doublet, truth.doublet)
    fp <- setdiff(pred.doublet, truth.doublet)
    tn <- intersect(pred.singlet, truth.singlet)
    fn <- setdiff(pred.singlet, truth.singlet)

    if(type=='precision'){
      precision <- length(tp)/(length(tp) + length(fp))
      acc.list[[method]] <- precision
    }
    if(type=='recall'){
      recall <- length(tp)/(length(tp) + length(fn))
      acc.list[[method]] <- recall
    }
    if(type=='TNR'){
      TNR <- length(tn)/(length(tn) + length(fp))
      acc.list[[method]] <- TNR
    }
  }
  return(acc.list)
}


#' Calculate precision, recall, or true negative rate (TNR) of identified doublets on multiple datasets.
#'
#' Calculate precision, recall, or TNR based on identified doublet indices and true doublet annotations on multiple datasets.
#' @param doublet.list.all A list of doublet indices for different doublet-detection methods on multiple datasets.
#' @param label.list A list of vectors of 0/1 doublet annotations on multiple datasets.
#'
#' @return A list of precision, recall, or TNR for different doublet-detection methods on multiple datasets.
#' @export
#'
#' @examples
#' precision.list.all <- FindACC.All(doublet.list.all = doublet.list.all, label.list = label.list, type = 'precision')
#' recall.list.all <- FindACC.All(doublet.list.all = doublet.list.all, label.list = label.list, type = 'recall')
#' tnr.list.all <- FindACC.All(doublet.list.all = doublet.list.all, label.list = label.list, type = 'TNR')
#'
FindACC.All <- function(doublet.list.all, label.list, type){
  acc.list.all <- list()
  for(i in 1:length(doublet.list.all)){
    doublet.list <- doublet.list.all[[i]]
    label <- label.list[[i]]
    data.name <- names(doublet.list.all)[i]; data.name
    cat('Calculate', type, 'on dataset:', data.name, '......\n')
    acc.list <- FindACC(doublet.list = doublet.list, label = label, type = type)
    acc.list.all[[data.name]] <- acc.list
  }
  return(acc.list.all)
}



#' Remove doublets for single method on single dataset
#'
#' Remove doublets based on identified doublet indices for single method on single dataset.
#' @param count A scRNA-seq data matrix.
#' @param label A vector of cell type annotations.
#' @param doublets A vector of identified doublet indices.
#'
#' @return A list of scRNA-seq data matrix and cell type annotations after removing doublets.
#' @export
#'
#' @examples
#'
RemoveDoublets <- function(count, label, doublets){
  count.removal <- count[,-doublets]; dim(count.removal); dim(count)
  label.removal <- label[-doublets]
  data.removal <- list(count=count.removal, label=label.removal)
  return(data.removal)
}


#' Remove doublets for multiple methods on single dataset
#'
#' Remove doublets based on identified doublet indices for multiple methods on single dataset.
#' @param count A scRNA-seq data matrix.
#' @param label A vector of cell type annotations.
#' @param doublet.list A list of identified doublet indices of different methods.
#'
#' @return A list of scRNA-seq data matrix and cell type annotations after removing doublets by different methods.
#' @export
#'
#' @examples
#' data.removal.list <- RemoveDoublets.Method(count = data.de$count, label = data.de$label.cluster, doublet.list = doublet.list)
#'
RemoveDoublets.Method <- function(count, label, doublet.list){
  data.removal.list <- list()
  for(method in names(doublet.list)){
    #method <- names(doublet.list)[1]; method
    data.removal <- RemoveDoublets(count = count, label = label, doublets = doublet.list[[method]])
    data.removal.list[[method]] <- data.removal
  }
  return(data.removal.list)
}


#' Remove doublets for multiple methods on multiple datasets
#'
#' Remove doublets based on identified doublet indices for different doublet-detection methods on multiple datasets.
#' @param count.list A list of scRNA-seq count matrices.
#' @param label.list A list doublet annoations (0/1).
#' @param doublet.list.all A list of identified doublet indices of different methods on multiple datasets.
#'
#' @return A list of scRNA-seq count matrices and doublet annotations after removing doublets.
#' @export
#'
#' @examples
#' data.removal.all <- RemoveDoublets.All(count.list=count.list, label.list=label.list, doublet.list.all=doublet.list.all)
#'
RemoveDoublets.all <- function(count.list, label.list, doublet.list.all){
  count.removal.list <- list()
  label.removal.list <- list()
  for(name in names(count.list)){
    #name <- names(count.list)[1]
    print(name)
    count <- count.list[[name]]; dim(count)
    label <- label.list[[name]]; table(label)
    for(method in names(doublet.list.all[[1]])){
      #method <- names(doublet.list.all[[1]])[1]
      print(method)
      doublets <- doublet.list.all[[name]][[method]]
      data.removal <- RemoveDoublets(count, label, doublets)
      count.removal.list[[name]][[method]] <- data.removal$count
      label.removal.list[[name]][[method]] <- data.removal$label
    }
  }
  data.removal.all <- list(count=count.removal.list, label=label.removal.list)
  return(data.removal.all)
}



#' Remove doublets for multiple methods on multiple datasets and doublet rates
#'
#' Remove doublets based on identified doublet indices for different doublet-detection methods on multiple datasets and doublet rates.
#' @param count.list A list of scRNA-seq count matrices.
#' @param label.list A list doublet annoations (0/1).
#' @param doublet.list.all.rate A list of identified doublet indices of different methods on multiple datasets and doublet rates.
#'
#' @return A list of scRNA-seq count matrices and doublet annotations after removing doublets based on multiple doublet rates.
#' @export
#'
#' @examples
#' doublet.list.all.rate <- FindDoublets.All.Rate(score.list.all = score.list.all, rates = seq(0.01, 0.25, 0.01))
#'
RemoveDoublets.All.Rate <- function(count.list, label.list, doublet.list.all.rate){
  data.removal.all.rate <- list()
  for(rate in names(doublet.list.all.rate)){
    print(rate)
    doublet.list.all <- doublet.list.all.rate[[rate]]
    data.removal.all <- RemoveDoublets.all(count.list, label.list, doublet.list.all)
    data.removal.all.rate[[rate]] <- data.removal.all
  }
  return(data.removal.all.rate)
}







