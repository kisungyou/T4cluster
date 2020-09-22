#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "utilities.h"
#include <cassert>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;
using namespace std;

// SECTION 1 : ELEMENTARY COMPUTATION ==========================================
// [[Rcpp::export]]
arma::mat cpp_shortestpath(arma::umat locs, arma::mat dists){
  // parameters
  unsigned int N = dists.n_rows;
  unsigned int K = dists.n_cols;
  // step 0. initialize
  arma::mat DIST(N,N); DIST.fill(arma::datum::inf);
  // step 1. assign with 'intersection' rule
  for (unsigned int n=0;n<N;n++){
    for (unsigned int k=0;k<K;k++){
      DIST(n,locs(n,k)-1) = dists(n,k);
    }
    DIST(n,n) = 0.0;
  }
  for (unsigned int i=0;i<(N-1);i++){
    for (unsigned int j=(i+1);j<N;j++){
      if (!(std::isfinite(DIST(i,j))&&std::isfinite(DIST(j,i)))){
        DIST(i,j) = arma::datum::inf;
        DIST(j,i) = arma::datum::inf;
      }
    }
  }
  // step 2. run iteration
  for (unsigned int k=0;k<N;k++){
    for (unsigned int i=0;i<N;i++){
      for (unsigned int j=0;j<N;j++){
        if (DIST(i,j) > DIST(i,k) + DIST(k,j)){
          DIST(i,j) = DIST(i,k) + DIST(k,j);
        }
      }
    }
  }
  return(DIST);
}
// [[Rcpp::export]]
arma::mat cpp_pdist(arma::mat X, int p){
  // prepare
  int N = X.n_rows;
  arma::mat output(N,N,fill::zeros);
  // iterate
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      output(i,j) = arma::norm(X.row(i)-X.row(j), p);
      output(j,i) = output(i,j);
    }
  }
  // return
  return(output);
}
// [[Rcpp::export]]
arma::mat cpp_pdist2(arma::mat X, arma::mat Y, int p){
  int M = X.n_rows;
  int N = Y.n_rows;
  
  arma::mat output(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      output(m,n) = arma::norm(X.row(m)-Y.row(n), p);
    }
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat cpp_pdistMP(arma::mat X, int p, int nCores){
  // prepare
  int N = X.n_rows;
  int d = X.n_cols;
  int useCores = 0;
  int detCores = 0;
  
  arma::mat output(N,N,fill::zeros);
#ifdef _OPENMP
  detCores = omp_get_num_procs();
  if (nCores <= 1){ // using 0 cores means I'll use 1 core
    useCores = 1;
  } else {
    if (detCores > nCores){
      useCores = nCores;
    } else {
      useCores = detCores;
    }
  }
#pragma omp parallel for num_threads(useCores) collapse(2) shared(output, N)
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      if (i<j){
        output(i,j) = arma::norm(X.row(i)-X.row(j),p);
        output(j,i) = output(i,j);  
      }
    }
  }
#else
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      output(i,j) = arma::norm(X.row(i)-X.row(j),p);
      output(j,i) = output(i,j);
    }
  }
#endif
  return(output);
}
// [[Rcpp::export]]
arma::uvec cpp_sample(int N, int m, arma::vec prob, bool replace){
  arma::uvec x     = arma::linspace<arma::uvec>(0L, N-1L, N);
  arma::vec myprob = prob/arma::accu(prob);
  arma::uvec output = Rcpp::RcppArmadillo::sample(x, m, replace, myprob);
  return(output);
}
// setdiff implementation
// https://stackoverflow.com/questions/29724083/trying-to-write-a-setdiff-function-using-rcpparmadillo-gives-compilation-error
// [[Rcpp::export]]
arma::uvec cpp_setdiff(arma::uvec& x, arma::uvec& y){
  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}


// SECTION 2 : K-MEANS AND GMM =================================================
// [[Rcpp::export]]
arma::urowvec label_kmeans(arma::mat data, int K, int maxiter){
  // parameters 
  int N = data.n_rows;
  
  // run k-means
  arma::mat means;
  bool status = arma::kmeans(means, arma::trans(data), K, random_subset, maxiter, false); // it returns (K x K) column means
  if (status == false){
    Rcpp::Rcout << "* k-means failed" << std::endl;
  }
  // need to compute pairwise distance matrix
  arma::mat kdist(K,N,fill::zeros);
  arma::colvec dcoli;
  for (int i=0; i<N; i++){
    dcoli = arma::trans(data.row(i));
    for (int j=0; j<K; j++){
      kdist(j,i) = arma::norm(means.col(j)-dcoli,2);
    }
  }
  urowvec gaus_ids = arma::index_min(kdist, 0);
  return(gaus_ids);
}
// [[Rcpp::export]]
arma::urowvec label_gmm(arma::mat data, int K, int maxiter){
  arma::gmm_full model;
  bool status = model.learn(data.t(), K, maha_dist, random_subset, maxiter, maxiter, 1e-10, false);
  if (status == false){
    Rcpp::Rcout << "* GMM failed" << std::endl;
  }
  urowvec gaus_ids = model.assign(data.t(), prob_dist);
  return(gaus_ids);
}

// SECTION 3 : BASIC SPECTRAL CLUSTERING =======================================
// [[Rcpp::export]]
Rcpp::List sc_unnormalized(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian
  arma::mat A = W; A.diag().fill(0.0);
  int  N = A.n_rows;
  arma::vec Dvec = arma::sum(A, 1);
  arma::mat Dmat = arma::diagmat(Dvec);
  arma::mat L = Dmat - A;
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, L);
  
  arma::mat dat = eigvec.head_cols(K); // (N x K) columns are smallest eigenvectors
  arma::urowvec output;
  if (usekmeans==true){
    output = label_kmeans(dat, K, maxiter);
  } else {
    output = label_gmm(dat, K, maxiter);
  }
  
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("embeds")=dat,
                            Rcpp::Named("labels")=output);
}
// [[Rcpp::export]]
Rcpp::List sc_normalNJW(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian
  arma::mat A = W; A.diag().fill(0.0);
  int N = A.n_rows;
  arma::vec Dvec = arma::sum(A, 1);
  arma::vec Dhalfinv(N,fill::zeros);
  double Dvalue = 0.0;
  for (int n=0;n<N;n++){
    Dvalue = Dvec(n);
    if (Dvalue > arma::datum::eps){
      Dhalfinv(n) = 1.0/std::sqrt(static_cast<float>(Dvalue));
    }
  }
  arma::mat Dhalfmat = arma::diagmat(Dhalfinv);
  arma::mat L = arma::eye(N,N) - Dhalfmat*A*Dhalfmat;
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, L);
  
  arma::mat dat = eigvec.head_cols(K); // (N x K) columns are smallest eigenvectors
  for (int n=0;n<N;n++){
    dat.row(n) = dat.row(n)/arma::norm(dat.row(n),2);
  }
  arma::urowvec output;
  if (usekmeans==true){
    output = label_kmeans(dat, K, maxiter);
  } else {
    output = label_gmm(dat, K, maxiter);
  }
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("embeds")=dat,
                            Rcpp::Named("labels")=output);
}
// [[Rcpp::export]]
Rcpp::List sc_normalSM(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian
  arma::mat A = W; A.diag().fill(0.0);
  int N = A.n_rows;
  arma::vec Dvec = arma::sum(A, 1);
  arma::vec Dinv(N,fill::zeros);
  double Dvalue = 0.0;
  for (int n=0;n<N;n++){
    Dvalue = Dvec(n);
    if (Dvalue > arma::datum::eps){
      Dinv(n) = 1.0/Dvalue;
    }
  }
  arma::mat Dinvmat = arma::diagmat(Dinv);
  arma::mat L = arma::eye(N,N) - Dinvmat*A;

  arma::cx_vec cxval;
  arma::cx_mat cxvec;
  arma::eig_gen(cxval, cxvec, L);
  arma::vec eigval = arma::real(cxval);
  arma::mat eigvec = arma::real(cxvec);
  
  arma::mat dat = eigvec.head_cols(K); // (N x K) columns are smallest eigenvectors
  arma::urowvec output;
  if (usekmeans==true){
    output = label_kmeans(dat, K, maxiter);
  } else {
    output = label_gmm(dat, K, maxiter);
  }
  
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("embeds")=dat,
                            Rcpp::Named("labels")=output);
}

// SECTION 4   : GMM-RELATED FUNCTIONS =========================================
// gmm_predict : prediction given new data and old model's fits 
//               be careful for transpose on X and Means
// [[Rcpp::export]]
arma::uvec gmm_predict(arma::mat X, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs){
  // model sizes
  int k = oldcovs.n_slices;
  int p = oldcovs.n_cols;
  
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldmeans)); // column centroids
  model.set_fcovs(oldcovs);
  model.set_hefts(arma::trans(oldweight));
  
  arma::uvec output = arma::trans(model.assign(arma::trans(X), prob_dist));
  return(output);
}
// gmm_sample    : sample from the gmm model 
// [[Rcpp::export]]
arma::mat gmm_sample(int n, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs){
  // model sizes
  int k = oldcovs.n_slices;
  int p = oldcovs.n_cols;
  
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldmeans)); // column centroids
  model.set_fcovs(oldcovs);
  model.set_hefts(arma::trans(oldweight));
  
  arma::mat output = arma::trans(model.generate(n));
  return(output);
}
// gmm_loglkd    : compute the log-likelihood of the data and the model ===
// [[Rcpp::export]]
double gmm_loglkd(arma::mat X, arma::colvec oldweight, arma::mat oldmeans, arma::cube oldcovs){
  // model sizes
  int k = oldcovs.n_slices;
  int p = oldcovs.n_cols;
  
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldmeans)); // column centroids
  model.set_fcovs(oldcovs);
  model.set_hefts(arma::trans(oldweight));
  
  double output = model.sum_log_p(arma::trans(X));
  return(output);
}