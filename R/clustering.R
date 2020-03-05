#' Selection of the number of cluster for the k-dissimilarities algorithm
#' 
#' It displays the value of the loss function / average silhouette width, for different values of \code{k}
#'
#'@importFrom Rcpp evalCpp sourceCpp
#'@importFrom cluster silhouette
#'@importFrom stats dist kmeans
#'@import ggplot2
#'@useDynLib GBClust
#'
#' @param D A n x n numeric matrix with the dissimilarities
#' @param k_max Maximum number of clusters to be considered
#' @param nstart Number of random initializations
#' @param method The graph that has to reported, i.e. \code{elbow} or \code{silhouette}
#' @return It return a \code{ggplot2} graph of the loss function / average silhouette width, for different values of \code{k}, from 1 up to \code{k_max}.
#' 
#' @export
#' 
kdiss_select <- function(D, k_max, nstart = 1, method="elbow"){
  
  n <- nrow(D)
  D <- matrix(D,n,n)
  loss <- numeric(k_max)
  sil  <- numeric(k_max-1)
  ncluster <- 2:k_max
  
  fit     <- kdiss(D = D, k = 1, nstart = nstart, trace=FALSE)
  loss[1] <- fit$loss
  
  for(k in ncluster){
    fit     <- kdiss(D = D, k = k, nstart = nstart, trace=FALSE)
    loss[k] <- fit$loss
    sil[k-1]  <- mean(silhouette(fit$cluster,D)[,3])
  }
  if(method=="elbow"){
    p <- ggplot(data=data.frame(ncluster=c(1,ncluster), loss = loss), aes(x=ncluster,y=loss)) + geom_point() + geom_line() + theme_bw() + xlab("Number of cluster") + ylab("Loss function")
  } 
  if(method=="silhouette"){
    p <- ggplot(data=data.frame(ncluster=ncluster, sil=sil), aes(x=ncluster,y=sil)) + geom_point() + geom_line() + theme_bw() + xlab("Number of cluster") + ylab("Average silhouette width")
  }
  p 
}

#' K-dissimilarities clustering
#' 
#' Perform the k-dissimilarities Algorithm described in Rigon, Herring and Dunson (2020).
#'
#' @param D A n x n numeric matrix with the dissimilarities
#' @param k The number of clusters to be considered 
#' @param nstart Number of random initializations
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced
#' @return 
#' \describe{
#'   \item{\code{cluster}}{The letters of the alphabet}
#'   \item{\code{loss}}{A vector of numbers}
#' }
#' 
#' @export
#' 
kdiss <- function(D, k, nstart = 1,  trace=FALSE){
  
  # Integrity checks
  n <- nrow(D)
  
  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)
  
  D <- matrix(D,n,n)
  
  if(k==1){
    G     <- rep(1,n)
    freq  <- as.numeric(table(G))
    best_fit  <- kdissimilarities_C(D = D, G = G,  freq = freq, trace = trace)
    return(best_fit)
  }
  
  # Random allocation with equal sizes
  G     <- sample(cut(seq(1,n),breaks=k,labels=FALSE))
  freq  <- as.numeric(table(G))
  best_fit  <- kdissimilarities_C(D = D, G = G,  freq = freq, trace = trace)
  if(nstart >= 2){
    for(r in 2:nstart){
      G     <- sample(cut(seq(1,n),breaks=k,labels=FALSE))
      freq  <- as.numeric(table(G))
      fit     <- kdissimilarities_C(D = D, G = G,  freq = freq, trace = trace)
      if(fit$loss < best_fit$loss) best_fit <- fit
    }
  }
  return(best_fit)
}

#' K-means clustering with uncertainty quantification
#' 
#' Perform the Gibbs-sampling for the k-means algorithm, as described in Rigon, Herring and Dunson. 
#'
#' @param x numeric matrix of the data
#' @param k The number of clusters to be considered. 
#' @param a_lambda Hyperparameter of the Gamma prior on the scale parameter
#' @param b_lambda Hyperparameter of the Gamma prior on on the scale parameter
#' @param R Number of MCMC samples after burn-in
#' @param burn_in Number of MCMC samples to be discarded as burn-in period
#' @param nstart Number of random initializations for the k-means algorithm
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return 
#' \describe{
#'   \item{\code{G}}{The letters of the alphabet}
#'   \item{\code{lambda}}{A vector of numbers}
#'   \item{\code{loss}}{A vector of numbers}
#'   \item{\code{G_map}}{A vector of numbers}
#'   \item{\code{loss_map}}{A vector of numbers}
#' }
#' 
#' @export
#' 
kmeans_gibbs <- function(x, k, a_lambda, b_lambda, R = 1000, burn_in = 1000, nstart=10, trace=FALSE){
  
  # Integrity checks
  n <- nrow(x)
  d <- ncol(x)
  
  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)
  
  if(trace){cat("Initialization of the algorithm\n")}
  
  fit_map  <- kmeans(x = x, centers = k, nstart = nstart, trace=FALSE)
  G_map    <- fit_map$cluster
  freq_map <- as.numeric(table(G_map))
  
  if(trace){cat("Starting the Gibbs sampling (R + burn-in) \n")}
  
  fit <- Gibbs_kmeans_C(R = R + burn_in, X = x, a_lambda = a_lambda, b_lambda = b_lambda, 
                      G = G_map, freq = freq_map, trace = trace)
  
  # Removing the burn-in
  fit$G      <- fit$G[-c(1:burn_in),]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss   <- fit$loss[-c(1:burn_in)]
  
  # Adding the MAP solution
  fit$G_map    <- G_map
  fit$loss_map <- fit_map$tot.withinss
  fit
}

#' K-dissimilarities clustering with uncertainty quantification
#' 
#' Perform the Gibbs-sampling for the k-dissimilarities algorithm
#'
#' @param x numeric matrix of of the data
#' @param k The number of clusters to be considered. 
#' @param p Power of the Minkowski distance
#' @param a_lambda Hyperparameter of the Gamma prior on the scale parameter
#' @param b_lambda Hyperparameter of the Gamma prior on on the scale parameter
#' @param R Number of MCMC samples after burn-in
#' @param burn_in Number of MCMC samples to be discarded as burn-in period
#' @param nstart Number of random initializations for the k-means algorithm
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return 
#' \describe{
#'   \item{\code{G}}{The letters of the alphabet}
#'   \item{\code{lambda}}{A vector of numbers}
#'   \item{\code{loss}}{A vector of numbers}
#'   \item{\code{G_map}}{A vector of numbers}
#'   \item{\code{loss_map}}{A vector of numbers}
#' }
#' 
#' @export
#' 
Minkowski_gibbs <- function(x, k, p, a_lambda, b_lambda, R = 1000, burn_in = 1000, nstart = 10,  trace=FALSE){
  
  # Integrity checks
  n <- nrow(x)
  d <- ncol(x)
  
  # Number of cluster must be smaller than n and greater or equal than 1
  stopifnot(k <= n)
  stopifnot(k >= 1)
  
  if(trace){cat("Initialization of the algorithm\n")}
  D <- as.matrix(dist(x, method="minkowski", p = p))

  fit_map <- kdiss(D = D, k = k, nstart = nstart, trace=FALSE)
  G_map    <- fit_map$cluster
  freq_map <- as.numeric(table(G_map))
  if(trace){cat("Starting the Gibbs sampling (R + burn-in) \n")}
  
  fit <- Gibbs_Mink_C(R = R + burn_in, D = D, d = d, a_lambda = a_lambda, b_lambda = b_lambda, 
                            G = G_map, freq = freq_map, trace = trace)
  
  # Removing the burn-in
  fit$G      <- fit$G[-c(1:burn_in),]
  fit$lambda <- fit$lambda[-c(1:burn_in)]
  fit$loss   <- fit$loss[-c(1:burn_in)]
  
  # Adding the MAP solution
  fit$G_map    <- G_map
  fit$loss_map <- fit_map$loss
  fit
}

#' Selection of the number of cluster for the k-dissimilarities algorithm
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
kmeans2_select <- function(x, k_max, nstart = 1, algorithm="kmeans") {
  
  n    <- nrow(x)
  p    <- ncol(x)
  x    <- matrix(x,n,p)
  loss <- numeric(k_max)
  ncluster <- 1:k_max
  
  for(k in ncluster){
    fit     <- kmeans2(x = x, k = k, nstart = nstart, algorithm=algorithm, trace=FALSE)
    loss[k] <- fit$loss
  }
  p <- ggplot(data=data.frame(ncluster=ncluster, loss = loss), aes(x=ncluster,y=loss)) + geom_point() + geom_line() + theme_bw() + xlab("Number of cluster") + ylab("Loss function")
  p 
}

#' K-Means^2 Clustering
#' 
#' Perform k-means and k-means^2 on a data matrix
#'
#' @param x numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param k The number of clusters to be considered. A random set of (distinct) rows in x is chosen as the initial centres.
#' @param nstart Number of random sets that has been chosen
#' @param algorithm The algorithm to be used
#' @param trace logical: if true, tracing information on the progress of the algorithm is produced.
#' @return 
#' \itemize{
#'   \item A - The letters of the alphabet.
#'   \item B - A vector of numbers.
#' }
#' 
#' @export
#' 
kmeans2 <- function(x, k, nstart = 1, algorithm="kmeans", trace=FALSE){
  
  # Integrity checks
  n <- nrow(x)
  
  stopifnot(k <= n)
  stopifnot(k >= 1)
  
  p <- ncol(x)
  x <- matrix(x,n,p)
  
  if(algorithm=="kmeans"){
    m         <- sample_centroids(x,k)
    best_fit  <- kmeans_C(x = x, m = m, k = k, trace = trace)
    if(nstart >= 2){
      for(r in 2:nstart){
        m         <- sample_centroids(x,k)
        fit     <- kmeans_C(x = x, m = m, k = k, trace = trace)
        if(fit$loss < best_fit$loss) best_fit <- fit
      }
    }
  }
  else if(algorithm=="kmeans2"){
    m         <- sample_centroids(x,k)
    best_fit  <- kmeans2_C(x = x, m = m, k = k, trace = trace)
    if(nstart >= 2){
      for(r in 2:nstart){
        m         <- sample_centroids(x,k)
        fit     <- kmeans_C(x = x, m = m, k = k, trace = trace)
        if(fit$loss < best_fit$loss) best_fit <- fit
      }
    }
  }
  
  return(best_fit)
}