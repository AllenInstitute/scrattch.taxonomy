% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{geneSymbolToEnsembl}
\alias{geneSymbolToEnsembl}
\title{Convert gene symbols to Ensembl IDs}
\usage{
geneSymbolToEnsembl(
  gene.symbols,
  ncbi.taxid = 9606,
  use.synonyms = TRUE,
  remove.duplicates = TRUE,
  includeNonMammalianSpecies = FALSE
)
}
\arguments{
\item{gene.symbols}{List of gene symbols to convert to Ensembl IDs.}

\item{ncbi.taxid}{The integer part of the NCBITaxon ID for the species want to convert genes between.}

\item{use.synonyms}{If TRUE (default) will search synonyms for current gene symbols to try and match Ensembl IDs}

\item{remove.duplicates}{If TRUE (default) any genes that share Ensembl IDs with any other genes will have their Ensembl IDs set to NA to avoid ambiguity. In cases where a duplicate is introduced in the synonyms, the original gene.symbol will retain the Ensembl IDs and synonym duplicates will be set to NA.}

\item{includeNonMammalianSpecies}{Default (FALSE) only considers mammalian species. Set to TRUE if non-mammalian species are considered (much slower).}
}
\value{
Ensembl IDs for the input gene.symbols (or NA if not found or if duplicated)
}
\description{
NOTE: This function requires an internet connection.
}
