#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* ROUTINES FOR FUNCTIONAL DATA ANALYSIS
 * (01) fpp_pdist_lp : Lp distance by vect (T), vecf (T x N)
 */

// (01) fpp_pdist_lp : Lp distance =============================================
double fpp_pdist_single(arma::vec xy, arma::vec tt, double p){
  int N = xy.n_elem;
  double output = 0.0;
  for (int n=0; n<(N-1); n++){
    output += (std::pow(xy(n),p)+std::pow(xy(n+1),p))*(tt(n+1)-tt(n))/2.0;
  }
  double outval = std::pow(output, 1.0/p);
  return(outval);
}

// [[Rcpp::export]]
arma::mat fpp_pdist_lp(arma::vec vect, arma::mat& vecf, double myp){
  // prepare
  int TT = vect.n_elem;
  int NN = vecf.n_cols;
  
  arma::vec vdiff(TT,fill::zeros);
  arma::mat output(NN,NN,fill::zeros);
  if (myp > 100){ // inf
    for (int i=0; i<(NN-1); i++){
      for (int j=(i+1); j<NN; j++){
        vdiff = arma::abs(vecf.col(i)-vecf.col(j));
        output(i,j) = vdiff.max();
        output(j,i) = output(i,j);
      }
    }
  } else {        // standard Lp
    for (int i=0; i<(NN-1); i++){
      for (int j=(i+1); j<NN; j++){
        vdiff = arma::abs(vecf.col(i)-vecf.col(j));
        output(i,j) = fpp_pdist_single(vdiff, vect, myp);
        output(j,i) = output(i,j);
      }
    }
  }
  return(output);
}