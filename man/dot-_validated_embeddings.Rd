% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkTaxonomy.R
\name{._validated_embeddings}
\alias{._validated_embeddings}
\title{This function will validate the embeddings in an AIT file}
\usage{
._validated_embeddings(AIT.anndata, messages, isValid, isWarning)
}
\arguments{
\item{AIT.anndata}{The AIT anndata object to validate.}

\item{messages}{The current messages to append to.}

\item{isValid}{The current isValid status.}

\item{isWarning}{The current isWarning status.}
}
\value{
A list with messages, an isValid logical call, and an isWarning logical call
}
\description{
This function will validate the embeddings in an AIT file
}
\keyword{internal}
