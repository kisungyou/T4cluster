#ifndef T4CLUSTER_UTILITIES_H
#define T4CLUSTER_UTILITIES_H

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// SECTION 1 : ELEMENTARY COMPUTATION
arma::mat  cpp_pdist(arma::mat X, int p);
arma::mat  cpp_pdist2(arma::mat X, arma::mat Y, int p);
arma::mat  cpp_pdistMP(arma::mat X, int p, int nCores);
arma::mat  cpp_shortestpath(arma::umat locs, arma::mat dists);
arma::uvec cpp_sample(int N, int m, arma::vec prob, bool replace);
arma::uvec cpp_setdiff(arma::uvec& x, arma::uvec& y);

// SECTION 2 : K-MEANS AND GMM
arma::urowvec label_kmeans(arma::mat data, int K, int maxiter);
arma::urowvec label_gmm(arma::mat data, int K, int maxiter);

// SECTION 3 : BASIC SPECTRAL CLUSTERING
// * use 'kmeans' or 'gmm' from armadillo with maxiter
Rcpp::List sc_unnormalized(arma::mat W, int K, bool usekmeans, int maxiter); // L=D-A
Rcpp::List sc_normalNJW(arma::mat W, int K, bool usekmeans, int maxiter);    // L=D^{-1/2}(D-A)D^{-1/2}
Rcpp::List sc_normalSM(arma::mat W, int K, bool usekmeans, int maxiter);     // L=D\(D-A)

// SECTION 4 : GMM-RELATED FUNCTIONS
arma::uvec gmm_predict(arma::mat X, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs);
arma::mat gmm_sample(int n, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs);
double gmm_loglkd(arma::mat X, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs);

// SECTION 5 : INTERNAL CRITERIA / CLUSTER VALIDITY INDEX
arma::mat cvi_helper_classmean(arma::mat X, arma::uvec label);   // compute class-wise mean
arma::field<arma::uvec> cvi_helper_classindex(arma::uvec label); // index for each label
int cvi_helper_nw(arma::uvec label);                             // number of pairs in the same cluster

#endif