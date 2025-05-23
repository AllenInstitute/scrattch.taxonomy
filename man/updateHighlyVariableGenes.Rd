% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{updateHighlyVariableGenes}
\alias{updateHighlyVariableGenes}
\title{Function to add or update highly variable genes for the current mode}
\usage{
updateHighlyVariableGenes(
  AIT.anndata,
  variable.genes = NULL,
  mode = AIT.anndata$uns$mode,
  list.name = paste0("highly_variable_genes_", AIT.anndata$uns$mode),
  default.list.name = "highly_variable_genes_standard"
)
}
\arguments{
\item{AIT.anndata}{A reference taxonomy anndata object.}

\item{variable.genes}{Set of variable genes to add to the list.name slot in obs. Can be provided as a logical or character vector.}

\item{list.name}{Which slot in obs to should the highly variable genes go?  Default is highly_variable_genes_\link{mode}}

\item{default.list.name}{Which slot should highly variable genes be copied from if none are provided? Default is highly_variable_genes_standard}
}
\value{
an AIT.anndata object with the updated/additional vector of highly variable genes
}
\description{
Function to add or update highly variable genes for the current mode
}
