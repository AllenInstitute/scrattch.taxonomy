% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{subsample_taxonomy}
\alias{subsample_taxonomy}
\title{Function to subsample cells from a taxonomy}
\usage{
subsample_taxonomy(cluster.names, cell_ids, dend = NULL, subsample = 2000)
}
\arguments{
\item{cluster.names}{A vector of cell type names}

\item{cell_ids}{A vector of cell ids to subsample from}

\item{dend}{A dendrogram object to use for subsampling}

\item{subsample}{The number of cells to retain per cluster (default = 2000)}
}
\value{
boolean vector for subsampling
}
\description{
Function to subsample cells from a taxonomy
}
\keyword{internal}
