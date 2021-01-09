/* 
 * METHODS FOR COMPUTING DISTANCES BETWEEN GAUSSIAN MIXTURES
 * (01) cpp_gmmdist_l2    : L2 distance
 * (02) cpp_gmmdist_base  : just compute pairwise distance using base metric
 * (03) cpp_gmmdist_cs    : cauchy-schwarz pdf divergence 
 * (04) cpp_gmmdist_jr    : jensen-renyi divergence
 * (05) cpp_gmmdist_tsl   : total square loss (bregman divergence)
 * (06) cpp_gmmdist_klga  : KL/Gaussian Approximation/Single Gaussian
 *      cpp_gmmdist_klsel : KL/Gaussian Approximation/Select
 * (07) cpp_gmmdist_klpog : KL/Product of Gaussians approximation
 * (08) cpp_gmmdist_klmb  : KL/Matched Bound Approximation
 * (09) cpp_gmmdist_klva  : KL/variational approximation
 * (10) cpp_gmmdist_klvub : KL/variational upper bound
 * (11) cpp_gmmdist_he    : Hilbert Embedding with Gaussian Kernel
 */


#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// weight   : length-k vector
// mean     : (k x p) matrix of row-stacked means
// variance : (p x p x k) cube for covariance matrices



// (01) cpp_gmmdist_l2 ---------------------------------------------------------
// [[Rcpp::export]]
double cpp_gmmdist_l2(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int N = weight1.n_elem;
  int M = weight2.n_elem;
  int p = mean1.n_cols; 
  
  // compute three terms
  arma::mat tmpcov(p,p,fill::zeros);
  double term1 = 0.0;
  double term2 = 0.0;
  double term3 = 0.0;
  
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      tmpcov = variance1.slice(i)+variance1.slice(j);
      term1 += 2.0*single_gaussian(mean1.row(i), mean1.row(j), tmpcov)*weight1(i)*weight1(j);
    }
  }
  for (int i=0; i<N; i++){
    tmpcov = variance1.slice(i)*2.0;
    term1 += single_gaussian(mean1.row(i), mean1.row(i), tmpcov)*weight1(i)*weight1(i);
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      tmpcov = variance2.slice(i)+variance2.slice(j);
      term2 += 2.0*single_gaussian(mean2.row(i), mean2.row(j), tmpcov)*weight2(i)*weight2(j);
    }
  }
  for (int i=0; i<M; i++){
    tmpcov = variance2.slice(i)*2.0;
    term2 += single_gaussian(mean2.row(i), mean2.row(i), tmpcov)*weight2(i)*weight2(i);
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<M; j++){
      tmpcov = variance1.slice(i) + variance2.slice(j);
      term3 += 2.0*single_gaussian(mean1.row(i), mean2.row(j), tmpcov)*weight1(i)*weight2(j);
    }
  }
  double output = std::sqrt(term1+term2-term3);
  return(output);
}

// (02) cpp_gmmdist_base : just compute pairwise distance using base metric ----
// [[Rcpp::export]]
arma::mat cpp_gmmdist_base(arma::mat mean1, arma::cube variance1, 
                           arma::mat mean2, arma::cube variance2, std::string basedist){
  int N = mean1.n_rows;
  int M = mean2.n_rows;
  int p = mean1.n_cols;
  
  arma::rowvec m1(p,fill::zeros);
  arma::rowvec m2(p,fill::zeros);
  arma::mat  s1(p,p,fill::zeros);
  arma::mat  s2(p,p,fill::zeros);
  arma::cube s2sqrt(p,p,M,fill::zeros);
  
  // special treatment : Wasserstein-2 geometry requires sqrtm(covariance)
  if (basedist=="wass2"){
    for (int m=0; m<M; m++){
      s2sqrt.slice(m) = arma::sqrtmat_sympd(variance2.slice(m));
    }
  }
  
  arma::mat output(N,M,fill::zeros);
  for (int i=0; i<N; i++){
    m1 = mean1.row(i);
    s1 = variance1.slice(i);
    for (int j=0; j<M; j++){
      m2 = mean2.row(j);
      s2 = variance2.slice(j);
      
      if (basedist=="l2"){
        output(i,j) = gauss2dist_l2(m1,s1,m2,s2);
      } else if (basedist=="wass2"){
        output(i,j) = gauss2dist_wass2(m1,s1,m2,s2,s2sqrt.slice(j));
      } else if (basedist=="cs"){
        output(i,j) = gauss2dist_cs(m1,s1,m2,s2);
      }
    }
  }
  
  return(output);
}

// (03) cpp_gmmdist_cs   : cauchy-schwarz pdf divergence -----------------------
// [[Rcpp::export]]
double cpp_gmmdist_cs(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int M = weight1.n_elem;
  int K = weight2.n_elem;
  int D = mean1.n_cols;
  double pi2term = std::pow(2.0*arma::datum::pi, static_cast<double>(D)/2.0);
  
  // computation : cross term
  double term1 = 0.0;
  for (int m=0; m<M; m++){
    for (int k=0; k<K; k++){
      term1 += weight1(m)*weight2(k)*single_gaussian(mean1.row(m),mean2.row(k),variance1.slice(m)+variance2.slice(k));
    }
  }
  // computation : first mixture terms
  double term2 = 0.0;
  for (int m=0; m<M; m++){
    term2 += std::pow(weight1(m), 2.0)/(pi2term*std::sqrt(arma::det(variance1.slice(m))));
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      term2 += 2.0*weight1(i)*weight1(j)*single_gaussian(mean1.row(i),mean1.row(j),(variance1.slice(i)+variance1.slice(j)));
    }
  }
  // computation : second mixture terms
  double term3 = 0.0;
  for (int k=0; k<K; k++){
    term3 += std::pow(weight2(k), 2.0)/(pi2term*std::sqrt(arma::det(variance2.slice(k))));
  }
  for (int i=0; i<(K-1); i++){
    for (int j=(i+1); j<K; j++){
      term3 += 2.0*weight2(i)*weight2(j)*single_gaussian(mean2.row(i),mean2.row(j),(variance2.slice(i)+variance2.slice(j)));
    }
  }
  
  // gather and return
  double output = -std::log(term1) + 0.5*std::log(term2) + 0.5*std::log(term3);
  return(output);
}

// (04) cpp_gmmdist_jr   : jensen-renyi divergence -----------------------------
double cpp_gmmdist_jrentropy(arma::vec pi, arma::mat mu, arma::cube sig){
  int N = sig.n_slices;
  double sumval = 0.0;
  for (int n=0; n<N; n++){ // diagonal/self terms
    sumval += pi(n)*pi(n)*single_gaussian(mu.row(n), mu.row(n), (2.0*sig.slice(n)));
  }
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      sumval += 2.0*pi(i)*pi(j)*single_gaussian(mu.row(i), mu.row(j), (sig.slice(i)+sig.slice(j)));
    }
  }
  double output = -std::log(sumval);
  return(output);
}
// [[Rcpp::export]]
double cpp_gmmdist_jr(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  arma::vec  weight3   = arma::join_vert(weight1, weight2)/2.0;
  arma::mat  mean3     = arma::join_vert(mean1, mean2);
  arma::cube variance3 = arma::join_slices(variance1, variance2);
 
 // now compute
 double Hfg = cpp_gmmdist_jrentropy(weight3, mean3, variance3);
 double Hf  = cpp_gmmdist_jrentropy(weight1, mean1, variance1);
 double Hg  = cpp_gmmdist_jrentropy(weight2, mean2, variance2);
 
 double output = Hfg - 0.5*(Hf+Hg);
 return(output);
}
// (05) cpp_gmmdist_tsl  : total square loss (bregman divergence) --------------
// [[Rcpp::export]]
double cpp_gmmdist_tsl(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                       arma::vec weight2, arma::mat mean2, arma::cube variance2){
  // parameters
  int M = weight1.n_elem;
  int N = weight2.n_elem;
  
  // compute : first term
  double d1 = 0.0;
  for (int m=0; m<M; m++){
    d1 += weight1(m)*weight1(m)*single_gaussian(mean1.row(m),mean1.row(m),(variance1.slice(m)*2.0));
  }
  for (int i=0; i<(M-1); i++){
    for (int j=(i+1); j<M; j++){
      d1 += 2.0*weight1(i)*weight1(j)*single_gaussian(mean1.row(i),mean1.row(j),(variance1.slice(i)+variance1.slice(j)));
    }
  }
  // compute : second term
  double d2 = 0.0;
  for (int n=0; n<N; n++){
    d2 += weight2(n)*weight2(n)*single_gaussian(mean2.row(n),mean2.row(n),(2.0*variance2.slice(n)));
  }
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      d2 += 2.0*weight2(i)*weight2(j)*single_gaussian(mean2.row(i),mean2.row(j),(variance2.slice(i)+variance2.slice(j)));
    }
  }
  // compute : cross term
  double d3 = 0.0;
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      d3 += weight1(m)*weight2(n)*single_gaussian(mean1.row(m),mean2.row(n),(variance1.slice(m)+variance2.slice(n)));
    }
  }
  // final output
  double term_top = d1 + d2 - 2.0*d3;
  double term_bot = std::sqrt(1.0 + 4.0*d2);
  double output   = term_top/term_bot;
  return(output);
}
// (06) cpp_gmmdist_klga : KL/Gaussian Approximation/Single Gaussian -----------
// [[Rcpp::export]]
double cpp_gmmdist_klga(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                        arma::vec weight2, arma::mat mean2, arma::cube variance2){
  int M = weight1.n_elem;
  int N = weight2.n_elem;
  int p = mean1.n_cols;
  
  arma::vec pi_f = weight1/arma::accu(weight1);
  arma::vec pi_g = weight2/arma::accu(weight2);
  
  arma::rowvec mu_diff(p,fill::zeros);
  arma::rowvec mu_f(p,fill::zeros);
  arma::rowvec mu_g(p,fill::zeros);
  arma::mat var_f(p,p,fill::zeros);
  arma::mat var_g(p,p,fill::zeros);
  
  // merge with 'f'
  for (int m=0; m<M; m++){
    mu_f = mu_f + pi_f(m)*mean1.row(m);
  }
  for (int m=0; m<M; m++){
    mu_diff = mu_f - mean1.row(m);
    var_f   = var_f + pi_f(m)*(variance1.slice(m) + arma::trans(mu_diff)*mu_diff);
  }
  // merge with 'g'
  for (int n=0; n<N; n++){
    mu_g = mu_g + pi_g(n)*mean2.row(n);
  }
  for (int n=0; n<N; n++){
    mu_diff = mu_g - mean2.row(n);
    var_g   = var_g + pi_g(n)*(variance2.slice(n) + arma::trans(mu_diff)*mu_diff);
  }
  
  // compute KL divergence
  double output = gauss2dist_kl(mu_f, var_f, mu_g, var_g);
  return(output);
}
// cpp_gmmdist_klsel : KL/Gaussian Approximation/Select
// [[Rcpp::export]]
double cpp_gmmdist_klsel(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                         arma::vec weight2, arma::mat mean2, arma::cube variance2){
  int M = weight1.n_elem;
  int N = weight2.n_elem;
  
  arma::mat distmat(M,N,fill::zeros);
  for (int m=0; m<M; m++){
    for (int n=0; n<N; n++){
      distmat(m,n) = gauss2dist_kl(mean1.row(m), variance1.slice(m), mean2.row(n), variance2.slice(n));
    }
  }
  
  double output = distmat.min();
  return(output);
}




// (11) cpp_gmmdist_he    : Hilbert Embedding with Gaussian Kernel -------------
// [[Rcpp::export]]
double cpp_gmmdist_he(arma::vec weight1, arma::mat mean1, arma::cube variance1,
                      arma::vec weight2, arma::mat mean2, arma::cube variance2,
                      double theta){
  // parameters
  int M = weight1.n_elem;
  int N = weight2.n_elem;
  int d = mean1.n_cols;
  double dd = static_cast<double>(d);
  double theta_d = std::pow(theta, dd);
  arma::mat theta_I = (theta*theta)*arma::eye<arma::mat>(d,d);
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat tmpcov(d,d,fill::zeros);
  arma::rowvec mudiff(d,fill::zeros);
  
  // iteration 1 : on the first component sets
  double term1 = 0.0;
  for (int m=0; m<M; m++){
    tmpcov = variance1.slice(m)*2.0 + theta_I;
    arma::eig_sym(eigval, eigvec, tmpcov);
    
  }
}
