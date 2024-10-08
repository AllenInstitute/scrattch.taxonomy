% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setupMapMyCells.R
\name{save_query_markers_to_uns}
\alias{save_query_markers_to_uns}
\title{Saves the contents of the query_markers.json file to AIT -> uns -> hierarchical.}
\usage{
save_query_markers_to_uns(AIT_anndata, query_markers_output_path)
}
\arguments{
\item{AIT_anndata}{AIT reference taxonomy object where the file contents will be saved.}

\item{query_markers_output_path}{Local file path to the generated or user provided query_markers.json file.}
}
\value{
AIT reference taxonomy with the query markers saved in uns -> hierarchical.
}
\description{
Saves the contents of the query_markers.json file to AIT -> uns -> hierarchical.
}
\keyword{internal}
