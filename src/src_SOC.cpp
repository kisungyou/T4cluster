#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/* SRC FUNCTIONS FOR LEARNING WITH CLUSTERINGS
 * (01) src_pcm : pairwise co-occurrence matrix
 * (02) src_psm : posterior similarity matrix
 */


// (01) src_pcm : pairwise co-occurrence matrix ================================
// single_coocurrence : given a row vector
arma::mat single_coocurrence(arma::irowvec partition){
  int N = partition.n_elem;
  arma::mat output(N,N,fill::eye);
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      if (partition(i)==partition(j)){
        output(i,j) = 1.0;
        output(j,i) = 1.0;
      }
    }
  }
  return(output);
}
// [[Rcpp::export]]
arma::mat src_pcm(arma::imat& clmat){
  int M = clmat.n_rows;
  int N = clmat.n_cols;
  
  arma::cube record(N,N,M,fill::zeros);
  for (int m=0;m<M;m++){
    record.slice(m) = single_coocurrence(clmat.row(m));
  }
  
  arma::mat output = arma::sum(record, 2);
  return(output);
}

// (02) src_psm : posterior similarity matrix ==================================
// [[Rcpp::export]]
arma::mat src_psm(arma::imat& clmat){
  int M = clmat.n_rows;
  int N = clmat.n_cols;
  
  arma::cube record(N,N,M,fill::zeros);
  for (int m=0;m<M;m++){
    record.slice(m) = single_coocurrence(clmat.row(m));
  }
  
  arma::mat output = arma::mean(record, 2);
  return(output);
}