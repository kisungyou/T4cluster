#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// compute W2 distance
double gauss_w2median_dist(arma::rowvec m1, arma::mat c1, arma::rowvec m2, arma::mat c2, arma::mat c2sqrt){
  arma::mat tmpmat = arma::sqrtmat_sympd(c2sqrt*c1*c2sqrt);
  double term1  = std::pow(arma::norm(m1-m2,2), 2.0);
  double term2  = arma::trace(c1 + c2 - 2.0*tmpmat);
  double output = std::sqrt(term1+term2);
  return(output);
}

// stack [mu, var] in horizontal way
arma::mat gauss_w2median_subproblem(arma::vec weight, arma::mat mean, arma::cube vars){
  // parameter & weight normalize
  int K = weight.n_elem;
  int P = mean.n_cols;
  arma::vec vecpi = weight/arma::accu(weight);
  
  // compute : mean
  arma::rowvec out_mean(P,fill::zeros);
  for (int k=0; k<K; k++){
    out_mean = out_mean + vecpi(k)*mean.row(k);
  }
  
  // compute : covariance
  // initialize
  arma::mat oldvar(P,P,fill::zeros);
  arma::mat newvar(P,P,fill::zeros);
  
  arma::mat oldvarsqrt(P,P,fill::zeros);
  arma::mat tmpvar(P,P,fill::zeros);
  for (int k=0; k<K; k++){
    oldvar = oldvar + vars.slice(k)*vecpi(k);
  }
  // iteration
  double incvar = 1000.0;
  int maxiter   = 200;
  double abstol = (1e-10)*static_cast<double>(P*P);
  for (int it=0; it<maxiter; it++){
    // compute sqrtmat for the old variance
    oldvarsqrt = arma::sqrtmat_sympd(oldvar);
    // fill newvar with zeros
    newvar.fill(0.0);
    // iterate over objects
    for (int k=0; k<K; k++){
      tmpvar = oldvarsqrt*vars.slice(k)*oldvarsqrt;
      newvar = newvar + (vecpi(k)*arma::sqrtmat_sympd(tmpvar));
    }
    // updater
    incvar = arma::norm(newvar-oldvar,"fro");
    oldvar = newvar;
    if (incvar < abstol){
      break;
    }
  }
  
  // return
  arma::mat output = arma::join_horiz(out_mean, oldvar);
  return(output);
}

// [[Rcpp::export]]
Rcpp::List gauss_w2median(arma::vec &weight, arma::mat &mean, arma::cube &vars, int maxiter=100, double abstol=1e-8){
  // parameter and weight normalization
  int K = weight.n_elem;
  int p = mean.n_cols;
  
  arma::vec pi_weight_old(K,fill::ones);
  arma::vec pi_weight_new(K,fill::ones);
  arma::vec pi_given = weight/arma::accu(weight);
  arma::vec pi_now(K,fill::zeros);
  arma::vec pi_new(K,fill::zeros);
  arma::mat tmp_run(p,p+1,fill::zeros);
  
  arma::rowvec now_mu(p,fill::zeros);
  arma::vec now_var(p,p,fill::zeros);
  
  // compute sqrtmat for variances
  arma::cube varsqrt(p,p,K,fill::zeros);
  for (int k=0; k<K; k++){
    varsqrt.slice(k) = arma::sqrtmat_sympd(vars.slice(k));
  }
  
  // iteration
  double inc_dist = 10000.0;
  for (int it=0; it<maxiter; it++){
    Rcpp::Rcout << "iteration " << it << " starts.." << std::endl;
    // update the current weight 
    for (int k=0; k<K; k++){
      pi_now(k) = pi_given(k)/pi_weight_old(k);
    }
    pi_now /= arma::accu(pi_now);
    
    // compute the new mean and variance
    tmp_run = gauss_w2median_subproblem(pi_now, mean, vars);
    now_mu  = arma::trans(tmp_run.col(0));
    now_var = tmp_run.tail_cols(p);
    
    // update the weight
    for (int k=0; k<K; k++){
      pi_weight_new(k) = gauss_w2median_dist(now_mu, now_var, mean.row(k), vars.slice(k), varsqrt.slice(k));
    }
    inc_dist = arma::norm(pi_weight_old-pi_weight_new,2);
    pi_weight_old = pi_weight_new;
    if (inc_dist < abstol){
      Rcpp::Rcout << "iteration " << it << " terminates.." << std::endl;
      break;
    }
  }
  
  // return the object
  Rcpp::List output;
  output["mean"] = now_mu;
  output["variance"] = now_var;
  return(output);
}