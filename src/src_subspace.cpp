#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::vec fast_loss_prj(int nS, int dS, int mS, arma::mat PS, arma::mat xS, arma::vec muS){
  int n = nS;
  double d  = static_cast<double>(dS);
  int m = mS;
  double m1 = static_cast<double>(n);
  double m2 = m1*m1;
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
