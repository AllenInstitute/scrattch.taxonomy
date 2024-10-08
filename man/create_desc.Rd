% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_utils.R
\name{create_desc}
\alias{create_desc}
\title{Create a generic description file}
\usage{
create_desc(
  dat,
  name = colnames(dat),
  use_label_columns = FALSE,
  start_columns = NULL
)
}
\arguments{
\item{dat}{any data frame that you would like to create a description file for}

\item{name}{desired names of each element in the description file (default is the column names)}

\item{use_label_columns}{should only columns containing "_label" be included (default = FALSE)}

\item{start_columns}{character vector of variables to include first in the list
(default = NULL, but "cluster" would be a common choice)}
}
\value{
a data.frame with columns "base", "name", and "type" for writing to tome
}
\description{
Create a generic description file
}
\keyword{internal}
