% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hyperparameter.R
\name{FindParameters}
\alias{FindParameters}
\title{Hyperparameter searching}
\usage{
FindParameters(
  count,
  label,
  method,
  type,
  n_neighbors = NULL,
  min_gene_variability_pctl = NULL,
  n_prin_comps = NULL,
  nfeatures = NULL,
  PCs = NULL,
  nf = NULL,
  includePCs = NULL,
  max_depth = NULL,
  k = NULL,
  d = NULL,
  ntop.cxds = NULL,
  ntop.bcds = NULL,
  n_components = NULL,
  n_top_var_genes = NULL,
  n_iters = NULL
)
}
\arguments{
\item{count}{A scRNA-seq count matrix.}

\item{label}{A vector of doublet annotations.}

\item{method}{A name of doublet-detection method.}

\item{type}{"AUPRC" or "AUROC".}

\item{n_neighbors}{Number of nearest neighbors in KNN classifier (Scrublet).}

\item{min_gene_variability_pctl}{The top percentile of highly variable genes (Scrublet).}

\item{n_prin_comps}{Number of principal components used to construct KNN classifer (Scrublet).}

\item{nfeatures}{Number of highly variable genes (DoubletFinder).}

\item{PCs}{Number of principal components used to construct KNN classifer (DoubletFinder).}

\item{nf}{Number of highly variable genes (scDblFinder).}

\item{includePCs}{Index of principal components to include in the predictors (scDblFinder).}

\item{max_depth}{Maximum depth of decision trees (scDblFinder).}

\item{k}{Number of nearest neighbors in KNN classifier (doubletCells).}

\item{d}{Number of principal components used to construct KNN classifer (doubletCells).}

\item{ntop.cxds}{Number of top variance genes to consider (cxds).}

\item{ntop.bcds}{Number of top variance genes to consider (bcds).}

\item{n_components}{Number of principal components used for clustering (DoubletDetection).}

\item{n_top_var_genes}{Number of highest variance genes to use (DoubletDetection).}

\item{n_iters}{Number of fit operations from which to collect p-values (DoubletDetection).}
}
\value{
A vector of hyperparameter combinations that maximize "AUPRC" or "AUROC".
}
\description{
Perform a grid search across user-specified hyperparameters.
Return a combination of hyperparameters that optimize AUPRC or AUROC.
}
\examples{
data.list <- ReadData(path = ".../real_datasets")
count.list <- data.list$count
label.list <- lapply(data.list$label, FUN = function(label){
 ifelse(label == 'doublet', 1, 0)
})
count <- count.list$`pbmc-1A-dm`
label <- label.list$`pbmc-1A-dm`
result.parameter.scDblFinder <- FindParameters(count, label, method = 'scDblFinder', type = 'AUPRC',
nf=c(500, 1000, 1500, 2000, 2500),
includePCs=c(3, 4, 5, 6, 7),
max_depth=c(3, 4, 5, 6, 7))

}
