% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KLdivergenceCallAndQC.R
\name{compute_KLdiv}
\alias{compute_KLdiv}
\title{Compute the KL Divergence}
\usage{
compute_KLdiv(
  query_probabilities,
  reference_probabilities,
  select.cl = NULL,
  select.cells = NULL
)
}
\arguments{
\item{query_probabilities}{A matrix of query cell mapping probabilities where rows represents cells, columns represent clusters and values rowsum to 1.}

\item{reference_probabilities}{A matrix of reference cluster probabilities where rows and columns both represents clusters values represent confusion between cluster mappings/clusterings and rowsum to 1.}

\item{select.cl}{An (optional) vector of cluster ids (e.g., row/col names of reference_probabilities) representing clusters for which KL divergence should be calculated using}

\item{select.cells}{An (optional) vector of cell ids (e.g., row names of query_probabilities) representing cells for which KL divergence should be calculated on}
}
\value{
A matrix of KL divergences for each requested cell (row) in each requested cluster (column)
}
\description{
From Wikipedia: "The Kullback–Leibler (KL) divergence (also called relative entropy and I-divergence) is a type of statistical distance that measures how one probability distribution P is different from a second, reference probability distribution Q. A simple interpretation of the KL divergence of P from Q is the expected excess surprise from using Q as a model when the actual distribution is P."  In scrattch.mapping this distance is required for computing categorical confidence calls used for assessing cell quality.
}
\keyword{external}
