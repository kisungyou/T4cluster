#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/* SRC FUNCTIONS FOR SUBSPACE CLUSTERING
 * (01) fast_loss_prj : compute loss function for MSM
 * (02) cpp_LRR       : LRR  simplest
 * (03) cpp_LRSC      : LRSC simplest
 * (04) cpp_EKSS_0    : zero iterations
 *      cpp_EKSS_T
 *      cpp_EKSS_affinity
 */

// (01) fast_loss_prj : compute loss function for MSM ==========================
// [[Rcpp::export]]
arma::vec fast_loss_prj(int nS, int dS, int mS, arma::mat PS, arma::mat xS, arma::vec muS){
  int n = nS;
  double d  = static_cast<double>(dS);
  // int m = mS;
  double m1 = static_cast<double>(n);
  // double m2 = m1*m1;
  arma::mat P  = PS;
  arma::mat x  = xS;
  arma::vec mu = muS;
  arma::vec kisung_vec;
  
  // compute distance
  arma::vec dist(n,fill::zeros);
  for (int i=0;i<n;i++){
    // dist(i) = d*std::pow(sum(m2*std::pow(P*(trans(x.row(i))-mu) - (trans(x.row(i))-mu),2)),0.5); // also very weird
    kisung_vec = P*(arma::trans(x.row(i))-mu) - (trans(x.row(i))-mu);
    dist(i) = d*m1*arma::norm(kisung_vec,2);
  }
  return(dist);
}

// (02) cpp_LRR : compute LRR ==================================================
//      input is directly given in (n x p) form
// [[Rcpp::export]]
Rcpp::List cpp_LRR(arma::mat& X, int par_k, int par_r){
  // compute
  // int N = X.n_rows;
  arma::mat U;
  arma::mat V;
  arma::vec s;
  arma::svd(U,s,V,X.t());
  
  // compute affinity
  arma::mat Vr = V.head_cols(par_r);
  arma::mat C  = Vr*Vr.t();
  arma::mat W  = arma::abs(C);
  W.diag().fill(0.0);
  // W /= W.max(); // this is something I put arbitrarily for making it as affinity
  
  // arma::vec Cvec(N,fill::zeros);
  // for (int n=0; n<N; n++){
  //   Cvec = arma::abs(C.col(n));
  //   C.col(n) = Cvec/arma::max(Cvec);
  // }
  // arma::mat W = (C + C.t())/2.0;
  // W.diag().fill(0.0);
  // 
  // run NJW spectral clustering & return
  bool my_kmeans = true;
  int my_maxiter = 100;
  Rcpp::List output = sc_normalNJW(W, par_k, my_kmeans, my_maxiter);
  return(output);
}

// (03) cpp_LRSC      : LRSC simplest ==========================================
//      input is directly given in (n x p) form
// [[Rcpp::export]]
Rcpp::List cpp_LRSC(arma::mat& X, int K, std::string algtype, double tau){
  // compute
  arma::mat U;
  arma::mat V;
  arma::vec s;
  arma::svd(U,s,V,X.t());
  
  // Rcpp::Rcout << "size of U : " << U.n_rows << "," << U.n_cols << std::endl;
  // Rcpp::Rcout << "size of V : " << V.n_rows << "," << V.n_cols << std::endl;
  // Rcpp::Rcout << "size of s : " << s.n_elem << std::endl;
  
  arma::mat C;
  // case branching
  if (algtype=="exact"){
    arma::mat V1 = V.cols(arma::find(s > arma::datum::eps));
    C = V1*V1.t();
  } else if (algtype=="relaxed"){
    double sqrtau = std::sqrt(tau);
    double sqrtauinv = 1.0/sqrtau;
    int N = s.n_elem;
    arma::vec s_trf(N,fill::zeros);
    for (int n=0; n<N; n++){
      if (s(n) > sqrtauinv){
        s_trf(n) = 1.0 - (1.0/(tau*s(n)*s(n)));
      }
    }
    C = V.head_cols(N)*arma::diagmat(s_trf)*arma::trans(V.head_cols(N));
  }
  
  // use the same approach as LRR
  arma::mat W  = arma::abs(C);
  W.diag().fill(0.0);
  
  bool my_kmeans = true;
  int my_maxiter = 100;
  Rcpp::List output = sc_normalNJW(W, K, my_kmeans, my_maxiter);
  return(output);
}


// (04) cpp_EKSS_0 & cpp_EKSS_T ================================================
arma::mat runif_stiefel(int p, int d){ // St(p,d) in R^{p x d}
  arma::mat X(p,d,fill::randn);
  arma::mat R = arma::inv_sympd(arma::sqrtmat_sympd(X.t()*X));
  return(X*R);
}
arma::uvec cpp_EKSS_label(arma::mat X, arma::cube projections){
  int N = X.n_rows;
  int d = projections.n_cols;
  int K = projections.n_slices;
  
  arma::mat scores(N,K,fill::zeros);
  arma::mat Xproj(N,d,fill::zeros);
  for (int k=0; k<K; k++){
    Xproj = X*projections.slice(k); // (N x d)
    for (int n=0; n<N; n++){
      scores(n,k) = arma::norm(Xproj.row(n), 2);
    }
  }
  
  arma::uvec label(N,fill::zeros);
  arma::rowvec score_vec(K,fill::zeros);
  for (int n=0; n<N; n++){
    score_vec = Xproj.row(n);
    label(n)  = score_vec.index_max();
  }
  return(label);
}
bool cpp_EKSS_not_K_vector(arma::uvec label, int K){
  arma::uvec uniques = arma::unique(label);
  arma::uvec idxnow;
  unsigned int Kuint = (unsigned int)K;
  if (uniques.n_elem < Kuint){
    return(true);
  } else {
    for (int k=0; k<K; k++){
      idxnow.reset();
      idxnow = arma::find(label==uniques(k));
      if (idxnow.n_elem < 2){
        return(true);
      }
    }
    return(false);
  }
}

// [[Rcpp::export]]
arma::uvec cpp_EKSS_0(arma::mat& X, int K, int d){
  // int N = X.n_rows;
  int p = X.n_cols;
  
  // generate random matrices from St(p,d)
  arma::cube Us(p,d,K,fill::zeros);
  for (int k=0; k<K; k++){
    Us.slice(k) = runif_stiefel(p,d);
  }
  
  // compute cluster label
  arma::uvec output = cpp_EKSS_label(X, Us);
  return(output);
}
// [[Rcpp::export]]
arma::uvec cpp_EKSS_T(arma::mat& X, int K, int d, int maxiter){
  int N = X.n_rows;
  int p = X.n_cols;
  
  // generate random matrices from St(p,d)
  arma::cube Us_old(p,d,K,fill::zeros);
  arma::cube Us_new(p,d,K,fill::zeros);
  for (int k=0; k<K; k++){
    Us_old.slice(k) = runif_stiefel(p,d);
  }
  
  arma::uvec label_old = cpp_EKSS_label(X, Us_old);
  if (cpp_EKSS_not_K_vector(label_old, K)){
    return(label_old);
  }
  arma::uvec label_new(N,fill::zeros);
  
  // main iteration
  arma::uvec vec_unique;
  arma::uvec idx_unique;
  arma::mat Xpart;
  arma::mat Xcov(p,p,fill::zeros);
  
  arma::vec eigval;
  arma::mat eigvec;
  for (int it=0; it<maxiter; it++){
    // find unique elements of the label vector
    vec_unique.reset();
    vec_unique = arma::unique(label_old);
    
    // run PCA for each cluster
    for (int k=0; k<K; k++){
      idx_unique.reset();
      Xpart.reset();
      eigval.reset();
      eigvec.reset();
      
      idx_unique = arma::find(label_old==vec_unique(k));
      Xpart = X.rows(idx_unique);
      Xcov  = arma::cov(Xpart);
      arma::eig_sym(eigval, eigvec, Xcov);
      
      Us_new.slice(k) = eigvec.tail_cols(d);
    }
    
    // get new label (+ conditioning)
    label_new = cpp_EKSS_label(X, Us_new);
    if (cpp_EKSS_not_K_vector(label_new, K)){
      return(label_old);
    } else {
      label_old = label_new;
      Us_old = Us_new;
    }
  }
  return(label_old);
}
// [[Rcpp::export]]
arma::mat cpp_EKSS_affinity(arma::umat &labels){ // these are column-stacked labels
  int N = labels.n_rows;
  int B = labels.n_cols;
  double BB = static_cast<double>(B);
  
  // construction of naive affinity matrix
  arma::mat A(N,N,fill::zeros);
  arma::uvec label_now(N,fill::zeros);
  for (int b=0; b<B; b++){
    label_now = labels.col(b);
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        if (label_now(i)==label_now(j)){
          A(i,j) = A(i,j) + (1.0/BB);
          A(j,i) = A(j,i) + (1.0/BB);
        }
      }
    }
  }
  return(A);
}