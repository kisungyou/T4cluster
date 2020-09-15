#ifndef T4CLUSTER_UTILITIES_H
#define T4CLUSTER_UTILITIES_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// SECTION 1 : ELEMENTARY COMPUTATION
arma::mat cpp_pdist(arma::mat X, int p);
arma::mat cpp_pdist2(arma::mat X, arma::mat Y, int p);
arma::mat cpp_pdistMP(arma::mat X, int p, int nCores);

// SECTION 2 : K-MEANS AND GMM
arma::urowvec label_kmeans(arma::mat data, int K, int maxiter);
arma::urowvec label_gmm(arma::mat data, int K, int maxiter);

// SECTION 3 : BASIC SPECTRAL CLUSTERING
// * use 'kmeans' or 'gmm' from armadillo with maxiter
Rcpp::List sc_unnormalized(arma::mat W, int K, bool usekmeans, int maxiter); // L=D-A
Rcpp::List sc_normalNJW(arma::mat W, int K, bool usekmeans, int maxiter);    // L=D^{-1/2}(D-A)D^{-1/2}
Rcpp::List sc_normalSM(arma::mat W, int K, bool usekmeans, int maxiter);     // L=D\(D-A)

#endif