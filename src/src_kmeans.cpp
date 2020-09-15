#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List arma_kmeans_random(arma::mat& X, int k, int maxiter){ // data should be stacked as columns
  // 1. compute means
  arma::mat means;
  bool status = arma::kmeans(means, X, k, random_subset, maxiter, false);
  if(status == false){
    Rcpp::stop("* alg.kmeans : Fitting k-means with random initialization failed.");
  }
  // 2. compute pairwise distance
  arma::mat pdmat = cpp_pdist2(arma::trans(X), arma::trans(means), 2);
  // 3. return
  return(Rcpp::List::create(Rcpp::Named("means")=arma::trans(means),
                            Rcpp::Named("pdmat")=pdmat));
}
// [[Rcpp::export]]
Rcpp::List arma_kmeans_kmeanspp(arma::mat& X, arma::mat &init, int k, int maxiter){ // data and initial mats are stacked as columns
  // parameter
  int N = X.n_cols;
  int P = X.n_rows;
  
  // copy
  arma::mat means(P,k,fill::zeros);
  for (int i=0;i<k;i++){
    means.col(i) = init.col(i);
  }
  // run
  bool status = arma::kmeans(means, X, k, keep_existing, maxiter, false);
  if(status == false){
    Rcpp::stop("* alg.kmeans : Fitting k-means with k-means++ initialization failed.");
  }
  // 2. compute pairwise distance
  arma::mat pdmat = cpp_pdist2(arma::trans(X), arma::trans(means), 2);
  // 3. return
  return(Rcpp::List::create(Rcpp::Named("means")=arma::trans(means),
                            Rcpp::Named("pdmat")=pdmat));
}