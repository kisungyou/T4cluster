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
arma::mat  gmm_sample(int n, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs);
double     gmm_loglkd(arma::mat X, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs);

// SECTION 5 : INTERNAL CRITERIA / CLUSTER VALIDITY INDEX
arma::mat cvi_helper_classmean(arma::mat X, arma::uvec label);   // compute class-wise mean
arma::field<arma::uvec> cvi_helper_classindex(arma::uvec label); // index for each label
int cvi_helper_nw(arma::uvec label);                             // number of pairs in the same cluster
int cvi_helper_nunique(arma::uvec label);               // number of unique elements in a label

// SECTION 6 : DISTANCE BETWEEN GAUSSIAN DISTRIBUTIONS
double single_gaussian(arma::rowvec x, arma::rowvec mu, arma::mat sig, bool logreturn=false);
double gauss2dist_l2(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // L2 Distance
double gauss2dist_wass2(arma::rowvec m1, arma::mat c1, arma::rowvec m2,              // 2-Wasserstein Distance
                        arma::mat c2, arma::mat c2sqrt);
double gauss2dist_cs(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // Cauchy-Schwarz Divergence
double gauss2dist_kl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // Kullback-Leibler Divergence
double gauss2dist_jr(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // Jensen-Renyi Divergence of Order 2
double gauss2dist_tsl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2); // total square loss (bregman divergence)
double gauss2dist_sl(arma::rowvec m1, arma::mat s1, arma::rowvec m2, arma::mat s2);  // square loss       (bregman divergence)

// SECTION 7 : GAUSSIAN DISTRIBUTION
arma::mat gauss_rmvnorm(int N, arma::vec mu, arma::mat var); // sample from a single gaussian

// SECTION 8 : NUMERICAL TOOLS
double integrate_1d(arma::vec &tseq, arma::vec &fval); // integrate 1-d discretized signal

#endif