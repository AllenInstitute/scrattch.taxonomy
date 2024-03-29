% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildTaxonomy.R
\name{buildTaxonomy}
\alias{buildTaxonomy}
\title{This function builds the minimum files required for Shiny}
\usage{
buildTaxonomy(
  counts,
  meta.data,
  feature.set,
  umap.coords,
  taxonomyDir = getwd(),
  taxonomyName = "AI_taxonomy",
  celltypeColumn = "cluster",
  cluster_colors = NULL,
  metadata_names = NULL,
  subsample = 2000,
  reorder.dendrogram = FALSE,
  dend = NULL
)
}
\arguments{
\item{counts}{A count matrix in sparse format: dgCMatrix.}

\item{meta.data}{Meta.data corresponding to count matrix. Rownames must be equal to colnames of counts. "clusters" must be provided (see celltypeColumn and notes).}

\item{feature.set}{Set of feature used to calculate dendrogram. Typically highly variable and/or marker genes.}

\item{umap.coords}{Dimensionality reduction coordiant data.frame with 2 columns. Rownames must be equal to colnames of counts.}

\item{taxonomyDir}{The location to save Shiny objects, e.g. "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/NHP_BG_20220104/"}

\item{taxonomyName}{The file name to assign for the Taxonomy h5ad.}

\item{celltypeColumn}{Column name corresponding to where the clusters are located (default="cluster")}

\item{cluster_colors}{An optional named character vector where the values correspond to colors and the names correspond to celltypes in celltypeColumn.  If this vector is incomplete, a warning is thrown and it is ignored. cluster_colors can also be provided in the metadata (see notes)}

\item{metadata_names}{An optional named character vector where the vector NAMES correspond to columns in the metadata matrix and the vector VALUES correspond to how these metadata should be displayed in Shiny. This is used for writing the desc.feather file later.}

\item{subsample}{The number of cells to retain per cluster}

\item{reorder.dendrogram}{Should dendogram attempt to match a preset order? (Default = FALSE).  If TRUE, the dendrogram attempts to match the celltype factor order as closely as possible (if celltype is a character vector rather than a factor, this will sort clusters alphabetically, which is not ideal).

Notes: Precomputed clusters must be provided.  In the anndata object these will be stored using the term "cluster".  If celltypeColumn is anything other than cluster, then any existing "cluster" column will be overwritten by celltypeColumn.  Values can be provided without colors and ids (e.g., "cluster") or with them (e.g., "cluster_label" + "cluster_color" + "cluster_id").  In this case cluster_colors is ignored and colors are taken directly from the metadata.  Cluster_id's will be overwritten to match dendrogram order.}

\item{dend}{Existing dendrogram associated with this taxonomy (e.g., one calculated elsewhere).  If NULL (default) a new dendrogram will be calculated based on the input \code{feature.set}}
}
\value{
AIT anndata object in the specified format (only if return.anndata=TRUE)
}
\description{
This function builds the minimum files required for Shiny
}
