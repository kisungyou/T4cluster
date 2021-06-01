#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* CORESET CONSTRUCTION
 * (1) coreset_18B : lightweight coreset 
 */

// (1) coreset_18B -------------------------------------------------------------
//     * draw size-M coreset
//     * apply k-means clustering with k  
// [[Rcpp::export]]
Rcpp::List coreset_18B(arma::mat& X, int K, int M, int maxiter){
  // PARAMETER
  int N = X.n_rows; double NN = static_cast<double>(N);

  // STEP 1 : COMPUTE MEAN AND DISTANCES
  arma::rowvec Xmean = arma::mean(X, 0);
  arma::vec    distsq(N,fill::zeros);
  double     dval = 0.0;
  for (int n=0; n<N; n++){
    dval = arma::norm(X.row(n)-Xmean,2);
    distsq(n) = dval*dval;
  }
  double distsqsum = arma::accu(distsq);
  
  // STEP 2 : COMPUTE PROBABILITY
  arma::vec probability(N,fill::zeros);
  for (int n=0; n<N; n++){
    probability(n) = (0.5/NN) + (0.5*distsq(n)/distsqsum);
  }
  
  // STEP 3 : DRAW INDEX FOR CORESET
  arma::uvec coreid = cpp_sample(N, M, probability, false);
  
  // STEP 4 : EXTRACT DATA AND SET WEIGHTS
  arma::mat core_data = X.rows(coreid);
  arma::vec core_weight = (1.0/probability(coreid))*(1.0/static_cast<double>(M));
  
  // RUN K-MEANS CLUSTERING
  arma::uvec sub_sample = cpp_sample(M,M,core_weight,true);
  arma::mat  sub_data = X.rows(sub_sample);
  arma::mat  kmeans;
  bool status = arma::kmeans(kmeans, arma::trans(sub_data), K, random_subset, maxiter, false);
  if(status == false){
    Rcpp::stop("* coreset18B routine failed.");
  }
  arma::mat  kcenters = arma::trans(kmeans);
  
  arma::mat  distmat(N,K,fill::zeros);
  for (int n=0; n<N; n++){
    for (int k=0; k<K; k++){
      distmat(n,k) = arma::norm(X.row(n)-kcenters.row(k), 2);
    }
  }
  arma::uvec cluster(N,fill::zeros);
  for (int n=0; n<N; n++){
    cluster(n) = arma::index_min(distmat.row(n));
  }
  
  // COMPUTE WCSS
  double wcss = 0.0;
  arma::vec distcol;
  arma::vec distsubvec;
  for (int k=0; k<K; k++){
    distcol    = distmat.col(k);
    distsubvec = distcol(arma::find(cluster==k));
    if (distsubvec.n_elem > 0){
      wcss += arma::dot(distsubvec, distsubvec);  
    }
  }
  
  // RETURN OUTPUT
  return(Rcpp::List::create(Rcpp::Named("means")=kcenters,
                            Rcpp::Named("cluster")=cluster,
                            Rcpp::Named("wcss")=wcss));
}


