#include "RcppArmadillo.h"
#include "utilities.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


/* CVI:Quality Index for the data-clustering match
 * (01) index_CH :
 */


// (01) index_CH ===============================================================
// [[Rcpp::export]]
double index_CH(arma::mat& X, arma::uvec label){
  // parameters
  arma::field<arma::uvec> classindex = cvi_helper_classindex(label);
  int N = X.n_rows;
  int K = classindex.n_elem;
  int P = X.n_cols;
  
  // prepare
  arma::uvec now_index;
  arma::rowvec vecdiff(P,fill::zeros);
  
  // compute class 
  arma::mat classmeans(K,P,fill::zeros);
  for (int k=0; k<K; k++){
    now_index.reset();
    now_index = classindex(k);
    if (now_index.n_elem < 2){
      classmeans.row(k) = X.row(now_index(0));
    } else {
      classmeans.row(k) = arma::mean(X.rows(now_index), 0);
    }
  }
  
  // compute Bg
  int now_size = 0;
  double Bg = 0.0;
  arma::rowvec datamean = arma::mean(X, 0);
  for (int k=0; k<K; k++){
    now_index.reset();
    now_index = classindex(k);
    vecdiff = classmeans.row(k)-datamean;
    Bg += static_cast<double>(now_index.n_elem)*arma::dot(vecdiff,vecdiff);
  }
  
  // compute Wg
  double Wg = 0.0;
  for (int k=0; k<K; k++){
    now_index.reset();
    now_index = classindex(k);
    now_size  = now_index.n_elem;
    for (int i=0; i<now_size; i++){
      vecdiff = classmeans.row(k) - X.row(now_index(i));
      Wg += arma::dot(vecdiff,vecdiff);
    }
  }
  
  // finalize
  double term1 = Bg/static_cast<double>(K-1);
  double term2 = Wg/static_cast<double>(N-K);
  return(term1/term2);
}
