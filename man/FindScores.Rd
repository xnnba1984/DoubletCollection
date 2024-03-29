% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindDoublets.R
\name{FindScores}
\alias{FindScores}
\title{Calculate doublet scores on single dataset}
\usage{
FindScores(
  count,
  methods,
  n_neighbors = round(0.5 * sqrt(dim(count)[2])),
  min_gene_variability_pctl = 85L,
  n_prin_comps = 30L,
  nfeatures = 2000,
  PCs = 10,
  nf = 1000,
  includePCs = 5,
  max_depth = 5,
  k = 50,
  d = 50,
  ntop.cxds = 500,
  ntop.bcds = 500,
  n_components = 30,
  n_top_var_genes = 10000,
  n_iters = 5
)
}
\arguments{
\item{count}{A scRNA-seq count matrix.}

\item{n_neighbors}{The number of nearest neighbors in KNN classifier (Scrublet).}

\item{min_gene_variability_pctl}{The top percentile of highly variable genes (Scrublet).}

\item{n_prin_comps}{Number of principal components used to construct KNN classifer (Scrublet).}

\item{nfeatures}{Number of highly variable genes (DoubletFinder).}

\item{PCs}{Number of principal components used to construct KNN classifer (DoubletFinder).}

\item{nf}{Number of highly variable genes (scDblFinder).}

\item{includePCs}{The index of principal components to include in the predictors (scDblFinder).}

\item{max_depth}{Maximum depth of decision trees (scDblFinder).}

\item{k}{The number of nearest neighbors in KNN classifier (doubletCells).}

\item{d}{Number of principal components used to construct KNN classifer (doubletCells).}

\item{ntop.cxds}{Number of top variance genes to consider (cxds).}

\item{ntop.bcds}{Number of top variance genes to consider (bcds).}

\item{n_components}{Number of principal components used for clustering (DoubletDetection).}

\item{n_top_var_genes}{Number of highest variance genes to use (DoubletDetection).}

\item{n_iters}{Number of fit operations from which to collect p-values (DoubletDetection).}

\item{method}{A name vector of doublet-detection methods.}
}
\value{
A list of doublet scores calculated by each doublet-detection method.
}
\description{
Call different computational doublet-detection methods to calculate doublet scores on single dataset.
}
\examples{
score.list <- FindScores(count = count.list$`J293t-dm`, methods = c('cxds','bcds','hybrid'))

}
