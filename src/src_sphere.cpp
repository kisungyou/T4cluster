#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// some geometric routines : rowvec format =====================================
// geometry_proj : Proj_x(z) = z - <x,z>x
// geometry_dist : acos(<x,y>)
// geometry_exp  : exponential map
// geometry_log  : logarithm map

arma::rowvec geometry_proj(arma::rowvec x, arma::rowvec z){
  arma::rowvec output = z - arma::dot(x,z)*x;
  return(output);
}
double geometry_dist(arma::rowvec x, arma::rowvec y){
  return(std::acos(static_cast<float>(arma::dot(x,y))));
}
arma::rowvec geometry_exp(arma::rowvec x, arma::rowvec u){
  double eps   = 0.00001;
  double unorm = arma::norm(u, 2);
  arma::rowvec output(x.n_elem, fill::zeros);
  if (unorm < eps){
    output = x;
  } else {
    output = std::cos(unorm)*x + (std::sin(unorm)/unorm)*u;
  }
  return(output);
}
arma::rowvec geometry_log(arma::rowvec x, arma::rowvec y){
  arma::rowvec projvec = geometry_proj(x, y-x);
  arma::rowvec output = (geometry_dist(x,y)/arma::norm(projvec, 2))*projvec;
  return(output);
}



// =============================================================================
/*
 * (1) sp_spkmeans  : basic algorithm for spherical k-means
 * (2) sp_gskmeans  : geodesic spherical k-means (mine!)
 * 
 * label : length-N vector of assignment 
 * index : length-K field of uvec's
 */

// (1) sp_spkmeans         : basic algorithm for spherical k-means =============
//     sp_spkmeans_centers : compute centers given index
//     sp_spkmeans_cost    : compute the objective
//     sp_spkmeans_label   : labeling update


arma::mat sp_spkmeans_centers(arma::mat X, arma::field<arma::uvec> index){
  // parameters
  unsigned int K = index.n_elem;
  unsigned int P = X.n_cols;
  
  // iterate
  arma::uvec cid;
  arma::mat output(K,P,fill::zeros);
  arma::rowvec tmpmean(P,fill::zeros);
  for (unsigned int k=0; k<K; k++){
    // clear & find
    cid.reset();
    tmpmean.fill(0.0);
    cid = index(k);
    
    // case branching
    if (cid.n_elem < 2){
      output.row(k) = X.row(cid(0));
    } else {
      tmpmean = arma::mean(X.rows(cid), 0);
      output.row(k) = tmpmean/arma::norm(tmpmean, 2);
    }
  }
  return(output);
}
double sp_spkmeans_cost(arma::mat X, arma::field<arma::uvec> index){
  // parameters
  unsigned int K = index.n_elem;
  
  // iterate
  double output = 0.0;
  for (unsigned int k=0; k<K; k++){
    output += arma::norm(arma::sum(X.rows(index(k)), 0),2);
  }
  
  return(output);
}
arma::uvec sp_spkmeans_label(arma::mat X, arma::mat centers){
  unsigned int N = X.n_rows;
  unsigned int K = centers.n_rows;
  
  arma::uvec output(N,fill::zeros);
  arma::vec rec_dist(K,fill::zeros);
  
  for (unsigned int n=0; n<N; n++){
    rec_dist.fill(0.0);
    for (unsigned int k=0; k<K; k++){
      rec_dist(k) = arma::dot(X.row(n), centers.row(k));
    }
    output(n) = arma::index_min(rec_dist);
  }
  return(output);
}


// [[Rcpp::export]]
Rcpp::List sp_spkmeans(arma::mat& X, int K, std::string initializer, int maxiter, double abstol, bool printer){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  
  // variables declaration
  arma::uvec old_label(N,fill::zeros);
  arma::uvec new_label(N,fill::zeros);
  arma::field<arma::uvec> old_index;
  arma::field<arma::uvec> new_index;
  arma::mat old_center(K,P,fill::zeros);
  arma::mat new_center(K,P,fill::zeros);
  double old_cost = 0.0;
  double new_cost = 10.0;
  double inc_cost = 100.0;
  
  // initialize : now we are just using GMM
  if (initializer=="gmm"){
    old_label = arma::trans(label_gmm(X, K, maxiter));
  } else if (initializer=="kmeans"){
    old_label = arma::trans(label_kmeans(X, K, maxiter));
  }
  
  old_label  = arma::trans(label_gmm(X, K, maxiter));
  if (cvi_helper_nunique(old_label) < K){ // no empty cluster allowed !
    Rcpp::stop("* spkmeans : initialization failed.");
  }
  old_index  = cvi_helper_classindex(old_label);
  old_center = sp_spkmeans_centers(X, old_index);
  old_cost   = sp_spkmeans_cost(X, old_index);

  // main iteration
  for (int iter=0; iter<maxiter; iter++){
    // Voronoi partitions
    new_label = sp_spkmeans_label(X, old_center);
    if (cvi_helper_nunique(new_label) < K){ // no empty cluster allowed !
      if (printer){
        Rcpp::Rcout << "* spkmeans : iteration " << iter+1 << " : clusters collapsed and stop now." << std::endl;
      }
      break;
    }
    new_index  = cvi_helper_classindex(new_label);
    new_center = sp_spkmeans_centers(X, new_index);
    new_cost   = sp_spkmeans_cost(X, new_index);
    
    // update
    inc_cost   = std::abs(new_cost-old_cost);
    old_cost   = new_cost;
    old_label  = new_label;
    old_index  = new_index;
    old_center = new_center;
    
    // cost change
    if (inc_cost < abstol){
      if (printer){
        Rcpp::Rcout << "* spkmeans : iteration " << iter+1 << " : convergence achieved." << std::endl;
      }
      break;
    }
    if (iter==(maxiter-1)){
      if (printer){
        Rcpp::Rcout << "* spkmeans : iteration " << iter+1 << " : convergence not achieved, maximum iteration reached." << std::endl;
      }  
    }
    if (printer){
      Rcpp::Rcout << "* spkmeans : iteration " << iter+1 << " : complete." << std::endl;
    }  
  }

  // return
  return(Rcpp::List::create(Rcpp::Named("means")=old_center,
                            Rcpp::Named("cluster")=old_label,
                            Rcpp::Named("cost")=old_cost));
}


// (2) sp_gskmeans  : geodesic spherical k-means (mine!) =======================
//     sp_gskmeans_frechet : compute Frechet mean
//     sp_gskmeans_centers : compute centers given index
//     sp_spkmeans_cost    : compute the objective
//     sp_spkmeans_label   : labeling update



arma::rowvec sp_gskmeans_frechet(arma::mat X){
  // parameters
  int maxiter = 100;
  double abstol = 0.0001;
  int P = X.n_cols;
  int N = X.n_rows;
  double dbN = static_cast<double>(N);
  
  // initialize
  arma::rowvec mean_old = arma::mean(X, 0);
  mean_old /= arma::norm(mean_old, 2);
  arma::rowvec mean_tmp(P,fill::zeros);
  arma::rowvec mean_new(P,fill::zeros);
  arma::rowvec grad_vec(P,fill::zeros);
  arma::mat logvecs(N,P,fill::zeros);
  double mean_inc = 1000.0;
  
  // iterate
  for (int it=0; it<maxiter; it++){
    // compute all the log-maps
    for (int n=0; n<N; n++){
      logvecs.row(n) = geometry_log(mean_old, X.row(n));
    }
    
    // compute the gradient
    grad_vec = arma::sum(logvecs, 0)/(dbN*2.0);
    
    // update
    mean_tmp = geometry_exp(mean_old, grad_vec);
    mean_new = mean_tmp/arma::norm(mean_tmp, 2);
    mean_inc = arma::norm(mean_new-mean_old, 2);
    mean_old = mean_new;
    
    // break
    if (mean_inc < abstol){
      break;
    }
  }
  
  // return
  return(mean_old);
}
arma::mat sp_gskmeans_centers(arma::mat X, arma::field<arma::uvec> index){
  // parameters
  unsigned int K = index.n_elem;
  unsigned int P = X.n_cols;
  
  // iterate
  arma::uvec cid;
  arma::mat output(K,P,fill::zeros);
  arma::rowvec tmpmean(P,fill::zeros);
  for (unsigned int k=0; k<K; k++){
    // clear & find
    cid.reset();
    tmpmean.fill(0.0);
    cid = index(k);
    
    // case branching
    if (cid.n_elem < 2){
      output.row(k) = X.row(cid(0));
    } else {
      output.row(k) = sp_gskmeans_frechet(X.rows(cid));
    }
  }
  return(output);
}
double sp_gskmeans_cost(arma::mat X, arma::mat centers, arma::field<arma::uvec> index){
  // parameters
  unsigned int K = index.n_elem;
  
  // iterate
  int Ksub = 0;
  arma::uvec cid;
  double output = 0.0;
  for (unsigned int k=0; k<K; k++){
    // clear & find
    cid.reset();
    cid = index(k);
    
    // sum it !
    if (cid.n_elem < 2){
      output += std::pow(geometry_dist(X.row(cid(0)), centers.row(k)), 2.0);
    } else {
      Ksub = cid.n_elem;
      for (int Kiter=0; Kiter<Ksub; Kiter++){
        output += std::pow(geometry_dist(X.row(cid(Kiter)), centers.row(k)), 2.0);
      }
    }
  }
  return(output);
}
arma::uvec sp_gskmeans_label(arma::mat X, arma::mat centers){
  unsigned int N = X.n_rows;
  unsigned int K = centers.n_rows;
  
  arma::uvec output(N,fill::zeros);
  arma::vec rec_dist(K,fill::zeros);
  
  for (unsigned int n=0; n<N; n++){
    rec_dist.fill(0.0);
    for (unsigned int k=0; k<K; k++){
      rec_dist(k) = geometry_dist(X.row(n), centers.row(k));
    }
    output(n) = arma::index_min(rec_dist);
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List sp_gskmeans(arma::mat& X, int K, std::string initializer, int maxiter, double abstol, bool printer){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  
  // variables declaration
  arma::uvec old_label(N,fill::zeros);
  arma::uvec new_label(N,fill::zeros);
  arma::field<arma::uvec> old_index;
  arma::field<arma::uvec> new_index;
  arma::mat old_center(K,P,fill::zeros);
  arma::mat new_center(K,P,fill::zeros);
  double old_cost = 0.0;
  double new_cost = 10.0;
  double inc_cost = 100.0;
  
  // initialize : now we are just using GMM
  if (initializer=="gmm"){
    old_label = arma::trans(label_gmm(X, K, maxiter));
  } else if (initializer=="kmeans"){
    old_label = arma::trans(label_kmeans(X, K, maxiter));
  }
  
  old_label  = arma::trans(label_gmm(X, K, maxiter));
  if (cvi_helper_nunique(old_label) < K){ // no empty cluster allowed !
    Rcpp::stop("* gskmeans : initialization failed.");
  }
  old_index  = cvi_helper_classindex(old_label);
  old_center = sp_gskmeans_centers(X, old_index);
  old_cost   = sp_gskmeans_cost(X, old_center, old_index);
  
  // main iteration
  for (int iter=0; iter<maxiter; iter++){
    // Voronoi partitions
    new_label = sp_gskmeans_label(X, old_center);
    if (cvi_helper_nunique(new_label) < K){ // no empty cluster allowed !
      if (printer){
        Rcpp::Rcout << "* gskmeans : iteration " << iter+1 << " : clusters collapsed and stop now." << std::endl;
      }
      break;
    }
    new_index  = cvi_helper_classindex(new_label);
    new_center = sp_gskmeans_centers(X, new_index);
    new_cost   = sp_gskmeans_cost(X, new_center, new_index);
    
    // update
    inc_cost   = std::abs(new_cost-old_cost);
    old_cost   = new_cost;
    old_label  = new_label;
    old_index  = new_index;
    old_center = new_center;
    
    // cost change
    if (inc_cost < abstol){
      if (printer){
        Rcpp::Rcout << "* gskmeans : iteration " << iter+1 << " : convergence achieved." << std::endl;
      }
      break;
    }
    if (iter==(maxiter-1)){
      if (printer){
        Rcpp::Rcout << "* gskmeans : iteration " << iter+1 << " : convergence not achieved, maximum iteration reached." << std::endl;
      }  
    }
    if (printer){
      Rcpp::Rcout << "* gskmeans : iteration " << iter+1 << " : complete." << std::endl;
    }  
  }
  
  // return
  return(Rcpp::List::create(Rcpp::Named("means")=old_center,
                            Rcpp::Named("cluster")=old_label,
                            Rcpp::Named("cost")=old_cost));
}