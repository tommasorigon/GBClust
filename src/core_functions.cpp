#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double loss_kmeans(const arma::mat& x_cluster){
  arma::mat x_bar = arma::mean(x_cluster);
  return(arma::accu(pow(x_cluster.each_row() - x_bar, 2)));
}

// [[Rcpp::export]]
double loss_average(const arma::mat& D_cluster){
  double out = 0.5*sum(arma::mean(D_cluster,1)); 
  return(out);
}

// [[Rcpp::export]]
double loss_kmeans2(const arma::mat& x_cluster){
  int nk = x_cluster.n_rows;
  arma::mat x_bar = arma::mean(x_cluster);
  return(nk*(log(arma::accu(pow(x_cluster.each_row() - x_bar, 2))) - log(nk)));
}


// [[Rcpp::export]]
List kdissimilarities_C(const arma::mat& D, arma::vec G, arma::vec freq, bool trace){
  
  // Initialization
  int n = D.n_rows;
  int K = freq.n_elem;
  
  // Output is a huge matrix
  arma::vec lprob(K);
  arma::uvec idx_cluster;
  arma::vec Losses(K);
  bool convergence = false;
  double total_loss_old = -1; 
  double total_loss;
  
  // Initialization of internal quantities
  bool skip;

  // Store the initial losses
  for(int j =0; j < K; j++){
    idx_cluster = find(G == j + 1); 
    Losses(j) = loss_average(D(idx_cluster,idx_cluster));
  }
  total_loss_old = sum(Losses);

  // Cycle of the Gibbs sampling
  while(!convergence){
    
    // Cycle of the observations
    for(int i = 0; i < n; i++) {
      freq(G(i)-1) = freq(G(i)-1) - 1;
      // Is the cluster of G_i a singleton?
      if(freq(G(i)-1) == 0){
        skip = true;
      } else {
        skip = false;
        // The ith element is not allocated
        idx_cluster = find(G == G(i)); // Note that here G(i) is equal to the previous value
        arma::vec diss = trans(D.row(i));
        Losses(G(i)-1) = (freq(G(i)-1) + 1)/freq(G(i)-1) *(Losses(G(i)-1) - mean(diss(idx_cluster)));
      }
      
      if(!skip){
        // Allocate the probability to move into another cluster
        
        G(i) = arma::datum::nan;
        arma::vec diss = trans(D.row(i));
        
        for(int j = 0; j < K; j++){
          // Identify the elements within the cluster
          idx_cluster = find(G == j + 1); // This  does not include G since it has been set to NaN
          lprob(j) =  freq(j)/(freq(j)+1)*mean(diss(idx_cluster)) - 1/(freq(j) + 1)*Losses(j);
        }
        // Allocate the value
        G(i) = index_min(lprob) + 1;
        
        // Number of distinct values remain the same, increase the frequency
        idx_cluster = find(G == G(i)); 
        Losses(G(i)-1) = freq(G(i)-1)/(freq(G(i)-1)+1) * Losses(G(i)-1) + mean(diss(idx_cluster));
      }
      
      freq(G(i)-1) = freq(G(i)-1) + 1;
    }
    
    // Store the initial losses
    total_loss=0;
    for(int j =0; j < K; j++){
      idx_cluster = find(G == j + 1); 
      Losses(j) = loss_average(D(idx_cluster,idx_cluster));
    }
    total_loss = sum(Losses);
    
    if(trace) {Rprintf("Loss function: %f \n", total_loss);}
  
    if((total_loss_old - total_loss) == 0){convergence=true;}
    
    total_loss_old = total_loss;
  }
  return(List::create(Named("cluster") = G, Named("loss") = total_loss));}

// [[Rcpp::export]]
List Gibbs_kmeans_C(int R, const arma::mat& X, double a_lambda, double b_lambda, arma::vec G, arma::vec freq, bool trace){
  
  // Initialization
  int n = X.n_rows;
  int d = X.n_cols;
  int K = freq.n_elem;
  double lambda;
  
  // How many times the output is displayed
  int R_show = floor(R/5);
  
  // Internal vectors and quantities
  arma::vec prob(K);  // Probabilities of allocation
  arma::vec lprob(K); // Log-probabilities of allocation
  IntegerVector clusters = Range(1,K); // Possible clusters
  bool skip;
  
  // Output is a huge matrix 
  arma::mat G_out(R,n);
  arma::vec lambda_out(R);
  arma::vec total_loss_out(R);
  arma::uvec idx_cluster;
  
  // Update the total_loss and the lambda parameter
  double total_loss=0;
  for(int j=0; j < K; j++){
    arma::uvec idx_cluster; idx_cluster = find(G == j + 1); // Now G_i has been
    total_loss = total_loss + loss_kmeans(X.rows(idx_cluster));
  }
  
  // Sample lambda from a Gamma distribution
  lambda = R::rgamma(a_lambda + 0.5*n*d, 1/(b_lambda + total_loss));
  
  // Cycle of the Gibbs Sampling
  for(int r = 0; r < R; r++){
    
    // Cycle of the Observations
    for(int i = 0; i < n; i++) {
      
      // Which is the old frequency?
      freq(G(i)-1) = freq(G(i)-1) - 1; // Reduces the associated frequencies
      
      // Is the cluster of G_i a singleton?
      if(freq(G(i)-1) == 0){
        skip = true;
      } else {
        skip = false;
      }
      
      if(!skip){
        // The ith element is not allocated
        G(i) = arma::datum::nan;
        // Compute the probability to move into another cluster
        for(int j = 0; j < K; j++){
          // Identify the elements within the cluster
          idx_cluster = find(G == j + 1); // This  does not include G since it has been set to NaN
          arma::uvec idx_cluster_i(freq(j)+1); idx_cluster_i.head(freq(j)) = idx_cluster; idx_cluster_i.tail(1) = i;
          lprob(j) = - lambda*(loss_kmeans(X.rows(idx_cluster_i)) - loss_kmeans(X.rows(idx_cluster)));
        }
        
        // Log-sum-exp trick
        lprob = lprob - max(lprob); prob  = exp(lprob); prob  = prob/sum(prob);
        // Sample the new value
        G(i) = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0];
      }
      // Update the frequency
      freq(G(i)-1) = freq(G(i)-1) + 1;
    }
    
    // Update the total_loss and the w parameter
    double total_loss=0;
    for(int j=0; j < K; j++){
      // Update the total loss
      arma::uvec idx_cluster; idx_cluster = find(G == j + 1); // Now G_i has been
      total_loss = total_loss + loss_kmeans(X.rows(idx_cluster));
    }
    
    // Sample w from a Gamma distribution
    lambda = R::rgamma(a_lambda + 0.5*n*d, 1/(b_lambda + total_loss));
    
    if(trace) {if((r+1)%R_show==0) {Rprintf("Iteration: %i \n", r+1);}}
    
    G_out.row(r) = trans(G);
    lambda_out(r) = lambda;
    total_loss_out(r) = total_loss;
  }
  return(List::create(Named("G") = G_out, Named("lambda") = lambda_out, Named("loss") = total_loss_out));
}

// [[Rcpp::export]]
List Gibbs_Mink_C(int R, const arma::mat& D, int d, double a_lambda, double b_lambda, arma::vec G, arma::vec freq, bool trace){
  
  // Initialization
  int n = D.n_rows;
  int K = freq.n_elem;
  double lambda;

  int R_show = floor(R/5);
  
  // Initialization of internal quantities
  bool skip;
  arma::vec Losses(K);
  arma::vec prob(K);  // Probabilities of allocation
  arma::vec lprob(K); // Log-probabilities of allocation
  IntegerVector clusters = Range(1,K); // Possible clusters
  
  // Output is a huge matrix 
  arma::mat G_out(R,n);
  arma::vec lambda_out(R);
  arma::vec total_loss_out(R);
  arma::uvec idx_cluster;

  // Store the initial losses
  
  for(int j =0; j < K; j++){
    idx_cluster = find(G == j + 1); 
    Losses(j) = loss_average(D(idx_cluster,idx_cluster));
  }
      
  double total_loss = sum(Losses);
  lambda = R::rgamma(a_lambda + n*d, 1/(b_lambda + total_loss));
  
  // Cycle of the Gibbs sampling
  for(int r = 0; r < R; r++){
    
  // Cycle of the observations
    for(int i = 0; i < n; i++) {
      // Which is the old frequency?
      
      // Which is the old frequency?
      freq(G(i)-1) = freq(G(i)-1) - 1; // Reduces the associated frequencies
      int freq_tmp = freq(G(i)-1);
      
      // Is the cluster of G_i a singleton?
      if(freq_tmp == 0){
        skip = true;
      } else if(freq_tmp > 0) {
        skip = false;
        // The ith element is not allocated
        arma::uvec idx_cluster; idx_cluster = find(G == G(i)); // Note that here G(i) is equal to the previous value
        arma::vec diss = trans(D.row(i));
        Losses(G(i)-1) = (freq(G(i)-1) + 1)/freq(G(i)-1) *(Losses(G(i)-1) - mean(diss(idx_cluster)));
      }

      if(!skip){
        // Allocate the probability to move into another cluster
        G(i) = arma::datum::nan;

        for(int j = 0; j < K; j++){
          // Identify the elements within the cluster
          idx_cluster = find(G == j + 1); // This  does not include G since it has been set to NaN
          arma::vec diss = trans(D.row(i));
          lprob(j) = - lambda*(freq(j)/(freq(j)+1)*mean(diss(idx_cluster)) - 1/(freq(j) + 1)*Losses(j));
        }
        // Log-sum-exp trick
        lprob = lprob - max(lprob); prob  = exp(lprob); prob  = prob/sum(prob);
        // Sample the new value
        G(i) = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0];
      
        // Number of distinct values remain the same, increase the frequency
        idx_cluster = find(G == G(i)); 
        arma::vec diss = trans(D.row(i));
        Losses(G(i)-1) = freq(G(i)-1)/(freq(G(i)-1)+1) * Losses(G(i)-1) + mean(diss(idx_cluster));
      }
      freq(G(i)-1) = freq(G(i)-1) + 1;
    }
    
    
    // Sample w from a Gamma distribution
    total_loss = sum(Losses);
    lambda = R::rgamma(a_lambda + n*d, 1/(b_lambda + total_loss));
    
    if(trace) {if((r+1)%R_show==0) {Rprintf("Iteration: %i \n", r+1);}}
    
    G_out.row(r) = trans(G);
    lambda_out(r) = lambda;
    total_loss_out(r) = total_loss;
  }
  return(List::create(Named("G") = G_out, 
                      Named("lambda") = lambda_out, 
                      Named("loss") = total_loss_out));
}


// [[Rcpp::export]]
arma::mat sample_centroids(const arma::mat& x, int k){
  
  
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat m(k,p); // Output
  
  // Initialization
  arma::uvec indexes(k);
  arma::vec prob(n, arma::fill::ones);
  IntegerVector clusters = Range(1,n); // Possible clusters
  
  double dist_min;
  double dist;
  
  // First centroid is selected randomly
  indexes(0) = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0] - 1;
  m.row(0) = x.row(indexes(0));
  
  // Sample the jth centroid
  for(int j=1; j<k; j++){
    // Compute the n probabilities
    for(int i=0; i<n; i++){
      // Compute the distance from the first centroid
      dist_min = sum(pow(x.row(i) - m.row(0), 2));
      // Compute the distance from the other centroids
      for(int h=1; h < j; h++){
        dist = sum(pow(x.row(i) - m.row(h), 2));
        if(dist < dist_min){dist_min = dist;}
      }
      // Store the probability
      prob(i) = dist_min;
    }
    // Sample the centroid
    prob = prob/sum(prob);
    indexes(j) = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0] - 1;
    m.row(j) = x.row(indexes(j));
  }
  return(m);
}

// [[Rcpp::export]]
List kmeans_C(const arma::mat& x, arma::mat m, int k, bool trace){
  
  // Initialization
  int n = x.n_rows;
  
  // Output is a huge matrix
  arma::vec G(n);
  arma::vec lprob(k);
  arma::uvec idx_cluster;
  arma::vec losses(k);
  bool convergence = false;
  double total_loss_old = -1; 
  double total_loss;
  
  // Iterations
  while(!convergence){
    
    // Cluster allocation
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < k; j++){
        losses(j) = sum(pow(x.row(i) - m.row(j), 2));
      }
      // Compute the probability to move into another cluster
      G(i) = index_min(losses) + 1;
    }
    
    total_loss = 0;
    for(int j = 0; j < k; j++){
      // Identify the elements within the cluster
      idx_cluster = find(G == j + 1); // This  does not include G since it has been set to NaN
      m.row(j) = mean(x.rows(idx_cluster));
      total_loss = total_loss + loss_kmeans(x.rows(idx_cluster));
    }
    
    if((total_loss_old - total_loss) == 0){convergence=true;} else {
      if(trace) {Rprintf("Loss function: %f \n", total_loss);}
      total_loss_old = total_loss;
    }
    
  }
  return(List::create(Named("cluster") = G, Named("centers") = m, Named("loss") = total_loss));
}

// [[Rcpp::export]]
List kmeans2_C(const arma::mat& x, arma::mat m, int k, bool trace){
  
  // Initialization
  int n = x.n_rows;
  int p = x.n_cols;
  
  // Output is a huge matrix
  arma::vec G(n);
  arma::vec sigma2(k, arma::fill::ones);
  arma::vec lprob(k);
  arma::mat x_cluster;
  arma::uvec idx_cluster;
  arma::vec losses(k);
  int nk;
  bool convergence = false;
  double total_loss_old = -1; 
  double total_loss;
  
  // Iterations
  while(!convergence){
    
    // Cluster allocation
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < k; j++){
        losses(j) = p*log(sigma2(j)) + sum(pow(x.row(i) - m.row(j), 2))/sigma2(j);
      }
      // Compute the probability to move into another cluster
      G(i) = index_min(losses) + 1;
    }
    
    total_loss = 0;
    for(int j = 0; j < k; j++){
      // Identify the elements within the cluster
      idx_cluster = find(G == j + 1); // This  does not include G since it has been set to NaN
      x_cluster = x.rows(idx_cluster);
      nk = x_cluster.n_rows;
      m.row(j)  = mean(x_cluster);
      sigma2(j) = arma::accu(pow(x_cluster.each_row() - m.row(j), 2))/(nk*p);
      total_loss = total_loss + loss_kmeans2(x_cluster);
    }
    
    if((total_loss_old - total_loss) == 0){convergence=true;} else {
      if(trace) {Rprintf("Loss function: %f \n", total_loss);}
      total_loss_old = total_loss;
    }
    
  }
  return(List::create(Named("cluster") = G, Named("centers") = m, Named("loss") = total_loss));
}
