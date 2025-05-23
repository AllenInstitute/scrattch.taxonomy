% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.checkBuildTaxonomyParams}
\alias{.checkBuildTaxonomyParams}
\title{Function to sanity check buildTaxonomy parameters}
\usage{
.checkBuildTaxonomyParams(
  counts,
  normalized.expr,
  meta.data,
  highly_variable_genes,
  marker_genes,
  embeddings,
  celltypeColumn,
  cluster_stats,
  taxonomyDir,
  title,
  dend
)
}
\arguments{
\item{counts}{A count matrix (cells x genes)}

\item{meta.data}{A data.frame with cell metadata}

\item{celltypeColumn}{The column name in meta.data that contains the cell type information}

\item{umap.coords}{A matrix of UMAP coordinates}
}
\value{
Stops the function if any of the parameters are not as expected
}
\description{
Function to sanity check buildTaxonomy parameters
}
\keyword{internal}
