% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{comp_medoids}
\alias{comp_medoids}
\title{Compute the medoids}
\usage{
comp_medoids(D, cluster)
}
\arguments{
\item{D}{A \code{n x n} numeric matrix containing the dissimilarities, i.e. the output of the functions \code{\link{dist}} or \code{\link{daisy}}.}

\item{cluster}{A clustering solution, i.e. the output of \code{\link{kdiss}}.}
}
\value{
\describe{
  \item{\code{medoids}}{Indexes of the medoids.}
}
}
\description{
Compute the medoids of a given clustering solution based on the corresponding dissimilarity matrix.
}
