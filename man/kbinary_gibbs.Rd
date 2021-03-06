% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{kbinary_gibbs}
\alias{kbinary_gibbs}
\title{K-binary Gibbs sampling}
\usage{
kbinary_gibbs(
  x,
  k,
  lambda = 1,
  R = 1000,
  burn_in = 1000,
  nstart = 10,
  trace = FALSE
)
}
\arguments{
\item{x}{binary matrix of the data.}

\item{k}{The number of clusters to be considered.}

\item{lambda}{Gibbs posterior tuning parameter.}

\item{R}{Number of MCMC samples after burn-in.}

\item{burn_in}{Number of MCMC samples to be discarded as burn-in period.}

\item{nstart}{Number of random initializations for the k-means algorithm.}

\item{trace}{logical: if true, tracing information on the progress of the algorithm is produced.}
}
\value{
\describe{
  \item{\code{G}}{A \code{R x n} matrix including the cluster labels for each MCMC iteration.}
  \item{\code{loss}}{A \code{R}-dimensional vector including the values of the loss function for each MCMC iteration.}
  \item{\code{G_map}}{Labels of the clusters at the lowest value of the posterior that has been computed.}
  \item{\code{loss_map}}{Lowest value of the loss that has been computed.}
}
}
\description{
Perform the Gibbs-sampling for the k-binary clustering. This function is complementary to \code{\link{kbinary}}, which is used instead to get a point estimate.
}
