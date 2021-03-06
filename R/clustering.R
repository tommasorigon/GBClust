#' Selection of the number of cluster for the k-dissimilarities algorithm
#'
#' It displays the value of the loss function / average silhouette width, for different values of \code{k}
#'
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom cluster silhouette
#' @importFrom stats dist kmeans
#' @import ggplot2
#' @useDynLib GBClust
#'
#' @param D A \code{n x n} numeric matrix containing the dissimilarities, typically the output of \code{\link{dist}} or \code{\link{daisy}}.
#' @param k_max Maximum number of clusters to be considered.
#' @param nstart Number of random initializations.
#' @param method The graph that will be displayed. Supported options are \code{method="elbow"}, which displays the loss function, or \code{method="silhouette"}. See \code{\link{silhouette}} for details about the latter.
#' @return It return a \code{\link{ggplot2}} graph of the loss function / average silhouette width, for \code{k=1,...,k_max}.
#'
#' @export
#'
kdiss_select <- function(D, k_max, nstart = 1, method = "elbow") {
  n <- nrow(D)
  D <- matrix(D, n, n)
  loss <- numeric(k_max)
  sil <- numeric(k_max - 1)
  ncluster <- 2:k_max

  fit <- kdiss(D = D, k = 1, nstart = nstart, trace = FALSE)
  loss[1] <- fit$loss

  for (k in ncluster) {
    fit <- kdiss(D = D, k = k, nstart = nstart, trace = FALSE)
    loss[k] <- fit$loss
    sil[k - 1] <- mean(silhouette(fit$cluster, D)[, 3])
  }
  if (method == "elbow") {
    p <- ggplot(data = data.frame(ncluster = c(1, ncluster), loss = loss), aes(x = ncluster, y = loss)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      xlab("Number of clusters") +
      ylab("Loss function")
  }
  if (method == "silhouette") {
    p <- ggplot(data = data.frame(ncluster = ncluster, sil = sil), aes(x = ncluster, y = sil)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      xlab("Number of clusters") +
      ylab("Average silhouette width")
  }
  p
}

#' K-dissimilarities algorithm
#'
#' Perform the so-called k-dissimilarities algorithm.
#'
#' @param D A \code{n x n} numeric matrix containing the dissimilarities, typically the output of \code{\link{dist}} or \code{\link{daisy}}.
#' @param k The number of clusters to be considered. See \code{\link{kdiss_select}} for selection criteria.
#' @param nstart Number of random initializations.
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return
#' \describe{
#'   \item{\code{cluster}}{Labels of the clusters at convergence}
#'   \item{\code{loss}}{Numeric value of the loss function at convergence}
#' }
#'
#' @export
#'
kdiss <- function(D, k, nstart = 1, trace = FALSE) {

  # Integrity checks
  n <- nrow(D)

  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)

  D <- matrix(D, n, n)

  if (k == 1) {
    G <- rep(1, n)
    freq <- as.numeric(table(G))
    best_fit <- kdissimilarities_C(D = D, G = G, freq = freq, trace = trace)
    return(best_fit)
  }

  # Random allocation with equal sizes
  G <- sample(cut(seq(1, n), breaks = k, labels = FALSE))
  freq <- as.numeric(table(G))
  best_fit <- kdissimilarities_C(D = D, G = G, freq = freq, trace = trace)
  if (nstart >= 2) {
    for (r in 2:nstart) {
      G <- sample(cut(seq(1, n), breaks = k, labels = FALSE))
      freq <- as.numeric(table(G))
      fit <- kdissimilarities_C(D = D, G = G, freq = freq, trace = trace)
      if (fit$loss < best_fit$loss) best_fit <- fit
    }
  }
  return(best_fit)
}

#' Compute the medoids
#'
#' Compute the medoids of a given clustering solution based on the corresponding dissimilarity matrix. 
#'
#' @param D A \code{n x n} numeric matrix containing the dissimilarities, i.e. the output of the functions \code{\link{dist}} or \code{\link{daisy}}.
#' @param cluster A clustering solution, i.e. the output of \code{\link{kdiss}}.
#' @return
#' \describe{
#'   \item{\code{medoids}}{Indexes of the medoids.}
#' }
#'
#' @export
#'
comp_medoids <- function(D, cluster) {

  # Integrity checks
  n <- nrow(D)
  k <- max(cluster)

  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)

  D <- matrix(D, n, n)
  medoids <- numeric(k)

  for (j in 1:k) {
    index <- which(cluster == j)
    n_k <- length(index)

    med <- index[1]
    loss <- sum(D[index[1], index])

    for (i in 2:n_k) {
      med_new <- index[i]
      loss_new <- sum(D[index[i], index])
      if (loss_new < loss) {
        med <- med_new
        loss <- loss_new
      }
    }
    medoids[j] <- med
  }
  medoids
}


#' K-dissimilarities (Minkowski) Gibbs sampling
#'
#' Perform the Gibbs-sampling for the k-dissimilarities clustering using the Minkowski distance. This function is complementary to \code{\link{kdiss}}, which is used instead to get a point estimate.
#'
#' @param x numeric matrix of of the data.
#' @param k The number of clusters to be considered.
#' @param p Power of the Minkowski distance.
#' @param a_lambda Hyperparameter of the Gamma prior on the scale parameter. The default \code{a_lambda = 0} leads to an improper prior.
#' @param b_lambda Hyperparameter of the Gamma prior on on the scale parameter. The default \code{a_lambda = 0} leads to an improper prior.
#' @param R Number of MCMC samples after burn-in.
#' @param burn_in Number of MCMC samples to be discarded as burn-in period.
#' @param nstart Number of random initializations for the \code{\link{kdiss}} algorithm, used to initialize the MCMC chain.
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return
#' \describe{
#'   \item{\code{G}}{Labels of the clusters at each \code{MCMC} iteration.}
#'   \item{\code{lambda}}{Numeric vector of the values of \code{lambda} at each MCMC iteration.}
#'   \item{\code{loss}}{Numeric vector of the loss function at each MCMC iteration.}
#'   \item{\code{G_map}}{Labels of the clusters obtained using \code{\link{kdiss}}, representing the maximum a posteriori.}
#'   \item{\code{loss_map}}{Numeric value of the loss function obtained using \code{\link{kdiss}}, representing the maximized loss.}
#' }
#' @export
#'
Minkowski_gibbs <- function(x, k, p, a_lambda = 0, b_lambda = 0, R = 1000, burn_in = 1000, nstart = 10, trace = FALSE) {

  # Integrity checks
  n <- nrow(x)
  d <- ncol(x)

  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)

  if (trace) {
    cat("Initialization of the algorithm\n")
  }
  D <- as.matrix(dist(x, method = "minkowski", p = p))

  fit_map <- kdiss(D = D, k = k, nstart = nstart, trace = FALSE)
  G_map <- fit_map$cluster
  freq_map <- as.numeric(table(G_map))
  if (trace) {
    cat("Starting the Gibbs sampling (R + burn-in) \n")
  }

  fit <- Gibbs_Mink_C(
    R = R + burn_in, D = D, d = d, a_lambda = a_lambda, b_lambda = b_lambda,
    G = G_map, freq = freq_map, trace = trace
  )

  # Removing the burn-in
  fit$G <- fit$G[-c(1:burn_in), ]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss <- fit$loss[-c(1:burn_in)]

  # Adding the MAP solution
  fit$G_map <- G_map
  fit$loss_map <- fit_map$loss
  fit
}

#' Selection of the number of cluster for the k-means algorithm
#'
#' It displays the value of the loss function for various choices of k
#'
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k_max The maximum number of clusters to be considered. A random set of (distinct) rows in x is chosen as the initial centres.
#' @param nstart Number of random sets that has been chosen
#' @param algorithm The algorithm to be used, either \code{kmeans} or \code{kmeans2}
#' @return It plots the loss function for different clustering solutions
#'
#' @export
#'
kmeans2_select <- function(x, k_max, nstart = 1, algorithm = "kmeans") {
  n <- nrow(x)
  p <- ncol(x)
  x <- matrix(x, n, p)
  loss <- numeric(k_max)
  ncluster <- 1:k_max

  for (k in ncluster) {
    fit <- kmeans2(x = x, k = k, nstart = nstart, algorithm = algorithm, trace = FALSE)
    loss[k] <- fit$loss
  }
  p <- ggplot(data = data.frame(ncluster = ncluster, loss = loss), aes(x = ncluster, y = loss)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    xlab("Number of cluster") +
    ylab("Loss function")
  p
}

#' K-means clustering
#'
#' Perform the k-means clustering on a data matrix.
#'
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k The number of clusters to be considered. A random set of (distinct) rows in x is chosen as the initial centers.
#' @param nstart Number of random sets that has been chosen.
#' @param algorithm The optimization algorithm to be used.
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return
#' \describe{
#'   \item{\code{cluster}}{Labels of the clusters at convergence.}
#'   \item{\code{centers}}{The value of the centroids at convergence.}
#'   \item{\code{loss}}{Numeric value of the loss function at convergence.}
#' }
#'
#' @export
#'
kmeans2 <- function(x, k, nstart = 1, algorithm = "kmeans", trace = FALSE) {

  # Integrity checks
  n <- nrow(x)

  stopifnot(k <= n)
  stopifnot(k >= 1)

  p <- ncol(x)
  x <- matrix(x, n, p)

  if (algorithm == "kmeans") {
    m <- sample_centroids(x, k)
    best_fit <- kmeans_C(x = x, m = m, k = k, trace = trace)
    if (nstart >= 2) {
      for (r in 2:nstart) {
        m <- sample_centroids(x, k)
        fit <- kmeans_C(x = x, m = m, k = k, trace = trace)
        if (fit$loss < best_fit$loss) best_fit <- fit
      }
    }
  }
  else if (algorithm == "kmeans2") {
    m <- sample_centroids(x, k)
    best_fit <- kmeans2_C(x = x, m = m, k = k, trace = trace)
    if (nstart >= 2) {
      for (r in 2:nstart) {
        m <- sample_centroids(x, k)
        fit <- kmeans_C(x = x, m = m, k = k, trace = trace)
        if (fit$loss < best_fit$loss) best_fit <- fit
      }
    }
  }

  return(best_fit)
}

#' K-means Gibbs sampling
#'
#' Perform the Gibbs-sampling for the k-means clustering. This function is complementary to \code{\link{kmeans2}}, which is used instead to get a point estimate.
#'
#' @param x A \code{n x d} numeric matrix of the data.
#' @param k The number of clusters to be considered.
#' @param a_lambda Hyperparameter of the Gamma prior on the scale parameter.
#' @param b_lambda Hyperparameter of the Gamma prior on on the scale parameter.
#' @param R Number of MCMC samples after burn-in.
#' @param burn_in Number of MCMC samples to be discarded as burn-in period.
#' @param nstart Number of random initializations for the k-means algorithm.
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return
#' \describe{
#'   \item{\code{G}}{A \code{R x n} matrix including the cluster labels for each MCMC iteration.}
#'   \item{\code{lambda}}{A \code{R}-dimensional vector including the values of lambda for each MCMC iteration.}
#'   \item{\code{loss}}{A \code{R}-dimensional vector including the values of the loss function for each MCMC iteration.}
#'   \item{\code{G_map}}{Labels of the clusters at the lowest value of the posterior that has been computed.}
#'   \item{\code{loss_map}}{Lowest value of the loss that has been computed.}
#' }
#'
#' @export
#'
kmeans_gibbs <- function(x, k, a_lambda, b_lambda, R = 1000, burn_in = 1000, nstart = 10, trace = FALSE) {

  # Integrity checks
  n <- nrow(x)
  d <- ncol(x)

  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)

  if (trace) {
    cat("Initialization of the algorithm\n")
  }

  fit_map <- kmeans(x = x, centers = k, nstart = nstart, trace = FALSE)
  G_map <- fit_map$cluster
  freq_map <- as.numeric(table(G_map))

  if (trace) {
    cat("Starting the Gibbs sampling (R + burn-in) \n")
  }

  fit <- Gibbs_kmeans_C(
    R = R + burn_in, X = x, a_lambda = a_lambda, b_lambda = b_lambda,
    G = G_map, freq = freq_map, trace = trace
  )

  # Removing the burn-in
  fit$G <- fit$G[-c(1:burn_in), ]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss <- fit$loss[-c(1:burn_in)]

  # Adding the MAP solution
  fit$G_map <- G_map
  fit$loss_map <- fit_map$tot.withinss
  fit
}


#' Selection of the number of cluster for the k-binary algorithm
#'
#' It displays the value of the loss function for various choices of k.
#'
#' @param x binary matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k_max The maximum number of clusters to be considered. A random set of (distinct) rows in x is chosen as the initial centers.
#' @param nstart Number of random sets that has been chosen.
#' @return It plots the loss function for different clustering solutions.
#'
#' @export
#'
kbinary_select <- function(x, k_max, nstart = 1) {
  n <- nrow(x)
  p <- ncol(x)
  x <- matrix(x, n, p)
  loss <- numeric(k_max)
  ncluster <- 1:k_max

  for (k in ncluster) {
    fit <- kbinary(x = x, k = k, nstart = nstart, trace = FALSE)
    loss[k] <- fit$loss
  }
  p <- ggplot(data = data.frame(ncluster = ncluster, loss = loss), aes(x = ncluster, y = loss)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    xlab("Number of clusters") +
    ylab("Loss function")
  p
}

#' K-binary clustering
#'
#' Perform the so-called k-binary clustering algorithm, for obtaining groups when the data are binary observations. 
#'
#' @param x binary matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k The number of clusters to be considered. A random set of (distinct) rows in x is chosen as the initial centers.
#' @param nstart Number of random sets that has been chosen.
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return
#' \describe{
#'   \item{\code{cluster}}{Labels of the clusters at convergence.}
#'   \item{\code{centers}}{The value of the centroids at convergence.}
#'   \item{\code{loss}}{Numeric value of the loss function at convergence.}
#' }
#'
#' @export
#'
kbinary <- function(x, k, nstart = 1, trace = FALSE) {

  # Integrity checks
  n <- nrow(x)

  stopifnot(k <= n)
  stopifnot(k >= 1)

  p <- ncol(x)
  x <- matrix(x, n, p)

  # Random allocation with equal sizes

  if (k == 1) {
    # Initialize randomly with Gibbs-sampling
    G <- rep(1, n)
    freq <- as.numeric(table(G))
    best_fit <- kbinary_C(x = x, k = k, G = G, freq = freq, trace = trace)
    return(best_fit)
  }

  # Random allocation with equal sizes
  G <- sample(cut(seq(1, n), breaks = k, labels = FALSE))
  freq <- as.numeric(table(G))

  best_fit <- kbinary_C(x = x, k = k, G = G, freq = freq, trace = trace)
  if (nstart >= 2) {
    for (r in 2:nstart) {
      G <- sample(cut(seq(1, n), breaks = k, labels = FALSE))
      freq <- as.numeric(table(G))
      fit <- kbinary_C(x = x, k = k, G = G, freq = freq, trace = trace)
      if (fit$loss < best_fit$loss) best_fit <- fit
    }
  }
  return(best_fit)
}

#' K-binary Gibbs sampling
#'
#' Perform the Gibbs-sampling for the k-binary clustering. This function is complementary to \code{\link{kbinary}}, which is used instead to get a point estimate.
#'
#' @param x binary matrix of the data.
#' @param k The number of clusters to be considered.
#' @param lambda Gibbs posterior tuning parameter.
#' @param R Number of MCMC samples after burn-in.
#' @param burn_in Number of MCMC samples to be discarded as burn-in period.
#' @param nstart Number of random initializations for the k-means algorithm.
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return
#' \describe{
#'   \item{\code{G}}{A \code{R x n} matrix including the cluster labels for each MCMC iteration.}
#'   \item{\code{loss}}{A \code{R}-dimensional vector including the values of the loss function for each MCMC iteration.}
#'   \item{\code{G_map}}{Labels of the clusters at the lowest value of the posterior that has been computed.}
#'   \item{\code{loss_map}}{Lowest value of the loss that has been computed.}
#' }
#'
#' @export
#'
kbinary_gibbs <- function(x, k, lambda = 1, R = 1000, burn_in = 1000, nstart = 10, trace = FALSE) {

  # Integrity checks
  n <- nrow(x)
  d <- ncol(x)

  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)

  if (trace) {
    cat("Initialization of the algorithm\n")
  }


  fit_map <- kbinary(x = x, k = k, nstart = nstart, trace = FALSE)
  G_map <- fit_map$cluster
  freq_map <- as.numeric(table(G_map))
  if (trace) {
    cat("Starting the Gibbs sampling (R + burn-in) \n")
  }

  fit <- Gibbs_kbinary_C(R = R + burn_in, X = x, G = G_map, freq = freq_map, lambda = lambda, trace = trace)

  # Removing the burn-in
  fit$G <- fit$G[-c(1:burn_in), ]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss <- fit$loss[-c(1:burn_in)]

  # Adding the MAP solution
  fit$G_map <- G_map
  fit$loss_map <- fit_map$loss
  fit
}
