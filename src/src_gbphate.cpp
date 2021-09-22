#include <RcppArmadillo.h>
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double loss_spkmeans(const arma::mat& x_cluster){
  arma::rowvec mean_euc = arma::mean(x_cluster, 0);
  arma::rowvec mean_sp  = mean_euc/arma::norm(mean_euc, 2);
  return(arma::accu(arma::pow(x_cluster.each_row() - mean_sp, 2)));
}

// MAIN FUNCTION on HYPERSPHERE ================================================
// [[Rcpp::export]]
Rcpp::List spkmeans_gibbs(int burn_in, int nsample, const arma::mat& X, double gamma_a, double gamma_b, arma::vec G, arma::vec freq, bool printer){
  // parameters
  int n = X.n_rows;
  int d = X.n_cols-1;
  int K = freq.n_elem;
  int R = burn_in + nsample;
  double lambda = 0.0;
  
  double nn = static_cast<double>(n);
  double dd = static_cast<double>(d);
  
  // internal vectors and quantities
  arma::vec prob(K);  // Probabilities of allocation
  arma::vec lprob(K); // Log-probabilities of allocation
  IntegerVector clusters = Range(1,K); // Possible clusters
  bool skip;
  
  // output is a huge matrix
  arma::mat G_out(R,n);
  arma::vec lambda_out(R);
  arma::vec total_loss_out(R);
  arma::uvec idx_cluster;
  
  // Update the total_loss and the lambda parameter
  double total_loss=0;
  for(int j=0; j < K; j++){
    idx_cluster.reset();
    idx_cluster = find(G == j + 1); // Now G_i has been
    total_loss = total_loss + loss_spkmeans(X.rows(idx_cluster));
  }
  
  // Sample lambda from a Gamma distribution
  lambda = R::rgamma(gamma_a + 0.5*nn*dd, 1/(gamma_b + total_loss));
  
  // Cycle of the Gibbs Sampling
  arma::uvec sampler;
  
  for (int r = 0; r < R; r++){
    
    // Cycle of the Observations
    for (int i = 0; i < n; i++) {
      
      // Which is the old frequency?
      freq(G(i)-1) = freq(G(i)-1) - 1; // Reduces the associated frequencies
      
      // Is the cluster of G_i a singleton?
      if(freq(G(i)-1) == 0){
        skip = true;
      } else {
        skip = false;
      }
      
      if(!skip){
        // The i-th element is not allocated
        G(i) = arma::datum::nan;
        // Compute the probability to move into another cluster
        for(int j = 0; j < K; j++){
          // Identify the elements within the cluster
          idx_cluster.reset();
          idx_cluster = find(G == j + 1); // This  does not include G since it has been set to NaN
          arma::uvec idx_cluster_i(freq(j)+1); 
          idx_cluster_i.head(freq(j)) = idx_cluster; 
          idx_cluster_i.tail(1) = i;
          lprob(j) = - lambda*(loss_spkmeans(X.rows(idx_cluster_i)) - loss_spkmeans(X.rows(idx_cluster)));
        }
        
        // Log-sum-exp trick
        lprob = lprob - max(lprob); 
        prob  = exp(lprob); 
        prob  = prob/sum(prob);
        
        // Sample the new value
        // G(i) = Rcpp::RcppArmadillo::sample(clusters, 1, TRUE, prob)[0];
        sampler.reset();
        sampler = cpp_sample(K, 1, prob, true);
        G(i) = sampler(0)+1;
      }
      // Update the frequency
      freq(G(i)-1) = freq(G(i)-1) + 1;
    }
    
    // Update the total_loss and the w parameter
    double total_loss=0;
    for(int j=0; j < K; j++){
      // Update the total loss
      idx_cluster.reset();
      idx_cluster = find(G == j + 1); // Now G_i has been
      total_loss = total_loss + loss_spkmeans(X.rows(idx_cluster));
    }
    
    // Sample w from a Gamma distribution
    lambda = R::rgamma(gamma_a + 0.5*nn*dd, 1/(gamma_b + total_loss));
    
    G_out.row(r) = trans(G);
    lambda_out(r) = lambda;
    total_loss_out(r) = total_loss;
    
    if (r%10 == 0){
      Rcpp::Rcout << "*         : Stage 2 - Iteration " << r+1 << "/" << R+1 << " complete." << std::endl;
    }
  }
  
  // return
  return(Rcpp::List::create(Named("G") = G_out.tail_rows(nsample), 
                            Named("lambda") = lambda_out.tail(nsample), 
                            Named("loss") = total_loss_out.tail(nsample)));
}

