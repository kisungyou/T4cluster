#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* GENERIC ROUTINES FOR CUSTOM GMM ROUTINES
 * eval_gaussian      : evaluate multivariate gaussian
 * eval_gaussian_single
 * eval_label         : my personalized labeling method
 * gmm_standard_gamma : update Gamma    (E-STEP)
 * gmm_standard_pi    : update Pi       (M-STEP)
 * gmm_standard_mean  : update Mean     (M-STEP)
 * gmm_standard_cov   : update Variance (M-STEP)
 * gmm_skeleton       : exemplary skeletal code for GMM
 * gmm_combine_wsum   : weighted sum
 * gmm_density        : density of a gmm model
 * gmm_pdist_wass2    : compute pairwise distance between GMM components (Wass2)
 * gmm_w2barycenter   : W2 barycenter of Gaussian distributions
 */

/* FUNCTIONS FOR GAUSSIAN MIXTURE MODELS
 * (01) gmm_armadillo : use Armadillo's native routine for Gaussian Mixture
 * (02) gmm_11R       : Regularized Covariance
 * (03) gmm_16Gfix    : use fixed weights
 * (04) gmm_03F       : single run of probability
 */
arma::vec eval_gaussian(arma::mat X, arma::rowvec mu, arma::mat Sig, bool logreturn=false){
  // parameters
  int n = X.n_rows; //double nn = static_cast<double>(n);
  int d = X.n_cols; double dd = static_cast<double>(d);
  
  // preparation
  double add1 = -(dd/2.0)*std::log(2.0*arma::datum::pi);
  double add2 = std::log(arma::det(Sig))*(-0.5);
  arma::vec outvec(n,fill::zeros);
  arma::rowvec xdiff(d,fill::zeros);
  arma::mat Sinv = arma::inv_sympd(Sig);
  for (int i=0;i<n;i++){
    xdiff = X.row(i) - mu;
    outvec(i) = -(arma::accu(xdiff*Sinv*xdiff.t())/2.0) + add1 + add2;
  }
  if (logreturn==true){
    return(outvec);
  } else {
    return(arma::exp(outvec));
  }
}
double eval_gaussian_single(arma::rowvec x, arma::rowvec mu, arma::mat sig, bool logreturn=false){
  int d = x.n_elem; double dd = static_cast<double>(d);
  double add1 = -(dd/2.0)*std::log(2.0*arma::datum::pi);
  double add2 = std::log(arma::det(sig))*(-0.5);
  arma::rowvec xdiff = x-mu;
  arma::mat Sinv = arma::pinv(sig);
  
  double output = -(arma::accu(xdiff*Sinv*xdiff.t())/2.0) + add1 + add2;
  if (logreturn==true){
    return(output);
  } else {
    return(std::exp(output));
  }
}

// [[Rcpp::export]]
arma::uvec eval_label(arma::mat& X, arma::mat parMU, arma::cube parSIG, arma::vec parPI){
  // parameters
  int N = X.n_rows;
  int K = parSIG.n_slices;
  // compute gamma
  arma::mat parGAMMA(N,K,fill::zeros);
  for (int k=0; k<K; k++){
    parGAMMA.col(k) = parPI(k)*eval_gaussian(X, parMU.row(k), parSIG.slice(k), false);
  }
  // for each component, find the maximal
  arma::uvec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::index_max(parGAMMA.row(n));
  }
  return(output);
}
arma::mat gmm_standard_gamma(arma::mat X, arma::mat kMu, arma::cube kSig, arma::vec kPi){
  // parameters
  int N = X.n_rows;
  int K = kPi.n_elem;
  int P = X.n_cols;
  
  arma::mat probmat(N,K,fill::zeros);
  arma::rowvec parmu(P,fill::zeros);
  arma::mat    parsig(P,P,fill::zeros);
  for (int k=0;k<K;k++){
    parmu  = kMu.row(k);
    parsig = kSig.slice(k);
    probmat.col(k) = eval_gaussian(X, parmu, parsig, false)*kPi(k);
  }
  
  arma::rowvec vecK(K,fill::zeros);
  for (int n=0;n<N;n++){
    vecK = probmat.row(n);
    probmat.row(n) = vecK/arma::accu(vecK);
  }
  return(probmat);
}
arma::vec gmm_standard_pi(arma::mat Gamma){
  // parameters
  int N = Gamma.n_rows; double NN = static_cast<double>(N);
  int K = Gamma.n_cols;
  arma::vec outPI(K,fill::zeros);
  for (int k=0;k<K;k++){
    outPI(k) = arma::accu(Gamma.col(k))/NN;
  }
  return(outPI);
}
arma::mat gmm_standard_mean(arma::mat X, arma::mat Gamma){
  // parameters
  int N = Gamma.n_rows;
  int K = Gamma.n_cols;
  int P = X.n_cols;
  
  // iterate
  arma::mat output(K,P,fill::zeros);
  arma::rowvec tmpvec(P,fill::zeros);
  for (int k=0;k<K;k++){
    tmpvec.fill(0.0);
    for (int n=0;n<N;n++){
      tmpvec += Gamma(n,k)*X.row(n);
    }
    output.row(k) = tmpvec/arma::accu(Gamma.col(k));
  }
  return(output);
}
arma::cube gmm_standard_cov(arma::mat X, arma::mat Gamma, arma::mat Mu, bool usediag=false){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  int K = Mu.n_rows;
  
  // iterate
  arma::cube output(P,P,K,fill::zeros);
  arma::mat  mslice(P,P,fill::zeros);
  arma::rowvec xdiff(P,fill::zeros);
  arma::mat  tmpmat(P,P,fill::zeros);
  double Nk = 0.0;
  for (int k=0;k<K;k++){
    // denominator
    Nk = arma::accu(Gamma.col(k));
    // numerator
    mslice.fill(0.0);
    for (int n=0;n<N;n++){
      xdiff   = X.row(n) - Mu.row(k);
      mslice += Gamma(n,k)*(xdiff.t()*xdiff);
    }
    tmpmat = mslice/Nk;
    if (usediag==true){
      output.slice(k) = arma::diagmat(tmpmat);
    } else {
      output.slice(k) = tmpmat;
    }
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List gmm_skeleton(arma::mat& X, int k){
  // PARAMETERS
  int n = X.n_rows; 
  int p = X.n_cols;
  int maxiter = 100;
  
  // PREPARATION : DEFINITION
  arma::mat oldGamma(n,k,fill::zeros);  // (NxK) class membership
  arma::mat newGamma(n,k,fill::zeros);
  arma::mat oldMU(k,p,fill::zeros);     // (KxP) row class means
  arma::mat newMU(k,p,fill::zeros);
  arma::cube oldSIG(p,p,k,fill::zeros); // (PxPxK) covariance slices
  arma::cube newSIG(p,p,k,fill::zeros);
  arma::vec oldPI(k,fill::zeros);       // (K) proportion
  arma::vec newPI(k,fill::zeros);
  double oldLKD = 0.0;
  double newLKD = 0.0;
  
  // PREPARATION : INITIALIZATION
  arma::urowvec initlabel = label_gmm(X, k, 5); // small number of simple task
  for (int i=0; i<n; i++){
    oldGamma(i,initlabel(i)) = 1.0;
  }
  oldPI  = gmm_standard_pi(oldGamma);
  oldMU  = gmm_standard_mean(X, oldGamma);
  oldSIG = gmm_standard_cov(X, oldGamma, oldMU, false); // full covariance
  oldLKD = gmm_loglkd(X, oldPI, oldMU, oldSIG);         // use armadillo routines
  
  // ITERATION
  for (int it=0; it<maxiter; it++){
    // {E} update gamma
    newGamma = gmm_standard_gamma(X, oldMU, oldSIG, oldPI);
    // {M} update parameters + compute loglkd
    newPI  = gmm_standard_pi(newGamma);
    newMU  = gmm_standard_mean(X, newGamma);
    newSIG = gmm_standard_cov(X, newGamma, newMU, false);
    newLKD = gmm_loglkd(X, newPI, newMU, newSIG);
    // Breaking Part
    if ((it>0)&&(newLKD <= oldLKD)){
      break;
    } else {
      oldGamma = newGamma;
      oldPI    = newPI;
      oldMU    = newMU;
      oldSIG   = newSIG;
      oldLKD   = newLKD; 
    }
  }
  
  // USE ARMADILLO ROUTINES
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldMU));
  model.set_fcovs(oldSIG);
  model.set_hefts(arma::trans(oldPI));
  
  // RETURN
  Rcpp::List output;
  output["means"]   = oldMU;
  output["covs"]    = oldSIG;
  output["weight"]  = oldPI;
  output["loglkd"]  = model.sum_log_p(arma::trans(X));
  output["cluster"] = arma::trans(model.assign(arma::trans(X), prob_dist));
  return(output);
}



// (01) gmm_armadillo ==========================================================
// [[Rcpp::export]]
Rcpp::List gmm_armadillo(arma::mat& X, int k, int maxiter, bool usediag){
  if (usediag==true){ // use diagonal
    arma::gmm_diag model;
    bool status = model.learn(arma::trans(X), k, maha_dist, random_subset, 10, maxiter, 1e-11, false);
    if (status==false){
      Rcpp::stop("* gmm : Fitting GMM with diagonal covariance failed.");
    } else {
      // successful, let's take out the elements and return
      
      arma::mat diagcovs = model.dcovs;
      int p = diagcovs.n_rows;
      int k = diagcovs.n_cols;
      arma::cube myfcovs(p,p,k,fill::zeros);
      for (int i=0; i<k; i++){
        myfcovs.slice(i) = arma::diagmat(diagcovs.col(i));
      }
      return Rcpp::List::create(Rcpp::Named("means")=arma::trans(model.means),
                                Rcpp::Named("covs")=myfcovs,
                                Rcpp::Named("weight")=model.hefts,
                                Rcpp::Named("loglkd")=model.sum_log_p(arma::trans(X)),
                                Rcpp::Named("cluster")=arma::trans(model.assign(arma::trans(X), prob_dist)));
    }
  } else {
    arma::gmm_full model;
    bool status = model.learn(arma::trans(X), k, maha_dist, random_subset, 10, maxiter, 1e-11, false);
    if (status==false){
      Rcpp::stop("* gmm : Fitting GMM with full covariance failed.");
    } else {
      // successful, let's take out the elements and return
      return Rcpp::List::create(Rcpp::Named("means")=arma::trans(model.means),
                                Rcpp::Named("covs")=model.fcovs,
                                Rcpp::Named("weight")=model.hefts,
                                Rcpp::Named("loglkd")=model.sum_log_p(arma::trans(X)),
                                Rcpp::Named("cluster")=arma::trans(model.assign(arma::trans(X), prob_dist)));
    }
  }
}

// (02) gmm_11R ================================================================
double gmm11R_objective(arma::mat S, arma::mat X, arma::mat Z, double lambda){
  double term1 = arma::trace(S*X);
  double term2 = log(arma::det(X));
  double term3 = lambda*arma::norm(vectorise(Z),1);
  double output = term1-term2+term3;
  return(output);
}
double gmm11R_shrinkage(double a, double kappa){
  double term1=0;
  double term2=0;
  if (a>kappa){    term1 = a-kappa;  }
  if (a<-kappa){   term2 = -a-kappa; }
  double output = term1-term2;
  return(output);
}
arma::mat gmm11R_ADMMprecision(arma::mat S, double lambda){
  // 1. parameters and set up
  const int max_iter  = 1000;
  const double abstol = 1e-6;
  const double reltol = 1e-3;
  const int    n      = S.n_cols;
  
  arma::mat X(n,n,fill::zeros);
  arma::mat X_hat(n,n,fill::zeros);
  arma::mat Z(n,n,fill::zeros);
  arma::mat Zold(n,n,fill::zeros);
  arma::mat U(n,n,fill::zeros);
  
  double rho   = 1.0;
  double alpha = 1.0;
  
  arma::colvec es(n,fill::zeros);
  arma::colvec xi(n,fill::zeros);
  arma::mat Q(n,n,fill::zeros);
  
  arma::vec objval(max_iter,fill::zeros);
  arma::vec r_norm(max_iter,fill::zeros);
  arma::vec s_norm(max_iter,fill::zeros);
  arma::vec eps_pri(max_iter,fill::zeros);
  arma::vec eps_dual(max_iter,fill::zeros);
  
  for (int k=0;k<max_iter;k++){
    // update X
    eig_sym(es,Q,rho*(Z-U)-S);
    for (int i=0;i<n;i++){
      xi(i) = (es(i)+sqrt(pow(es(i),2)+4*rho))/(2*rho);
    }
    X = Q*arma::diagmat(xi)*Q.t();
    
    // update Z with relaxation
    Zold = Z;
    X_hat = alpha*X + (1-alpha)*Zold;
    for (int i=0;i<n;i++){
      for (int j=0;j<n;j++){
        Z(i,j) = gmm11R_shrinkage(X_hat(i,j)+U(i,j), lambda/rho);
      }
    }
    
    // update U
    U = U + (X_hat-Z);
    
    // diagnostics
    objval(k) = gmm11R_objective(S, X, Z, lambda);
    r_norm(k) = arma::norm(X-Z,"fro");
    s_norm(k) = arma::norm(-rho*(Z-Zold),"fro");
    if (norm(X,"fro")>norm(Z,"fro")){
      eps_pri(k) = static_cast<double>(n)*abstol + reltol*norm(X,"fro");
    } else {
      eps_pri(k) = static_cast<double>(n)*abstol + reltol*norm(Z,"fro");
    }
    eps_dual(k) = static_cast<double>(n)*abstol + reltol*norm(rho*U, "fro");
    
    if ((r_norm(k)<eps_pri(k))&&(s_norm(k)<eps_dual(k))){
      break;
    }
  }
  return(X);
}
arma::cube gmm11R_precision(arma::mat X, arma::mat Gamma, arma::mat Mu, double lambda){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  int K = Mu.n_rows;
  
  arma::mat  Ak(P,P,fill::zeros);
  arma::rowvec xdiff(P,fill::zeros);
  arma::rowvec muk(P,fill::zeros);
  double Nk = 0.0;
  arma::cube output(P,P,K,fill::zeros);
  for (int k=0;k<K;k++){
    // initialize
    Ak.fill(0.0);
    Nk  = arma::accu(Gamma.col(k));
    muk = Mu.row(k);
    for (int n=0;n<N;n++){
      xdiff = X.row(n) - muk;
      Ak   += (Gamma(n,k)/Nk)*(xdiff.t()*xdiff);
    }
    // compute
    output.slice(k) = gmm11R_ADMMprecision(Ak, lambda); // precision to covariance
  }
  return(output);
}
// [[Rcpp::export]]
Rcpp::List gmm_11R(arma::mat& X, int K, double lambda, int maxiter, bool usediag){
  // parameters
  int N = X.n_rows; // double NN = static_cast<double>(N);
  int P = X.n_cols;
  int mit = maxiter; // max iter
  
  // preparation : definition
  arma::mat oldGamma(N,K,fill::zeros);  // (NxK) class membership
  arma::mat newGamma(N,K,fill::zeros);
  arma::mat oldMU(K,P,fill::zeros);     // (KxP) row class means
  arma::mat newMU(K,P,fill::zeros);
  arma::cube oldSIG(P,P,K,fill::zeros); // (PxPxK) covariance slices
  arma::cube newSIG(P,P,K,fill::zeros);
  arma::cube oldPRE(P,P,K,fill::zeros);
  arma::cube newPRE(P,P,K,fill::zeros);
  arma::vec oldPI(K,fill::zeros);       // (K) proportion
  arma::vec newPI(K,fill::zeros);
  double oldLKD = 0.0;
  double newLKD = 0.0;
  
  int par0 = (P*K) + (K-1); // mean + proportion
  int parS = 0;             // covariance part
  int oldPar = 0;
  int newPar = 0;
  
  double valijk = 0.0;
  double thres0 = arma::datum::eps;
  
  // preparation : initialize
  arma::urowvec initlabel = label_gmm(X, K, 5); // small number of simple task
  for (int n=0; n<N; n++){
    oldGamma(n,initlabel(n)) = 1.0;
  }
  oldPI  = gmm_standard_pi(oldGamma);
  oldMU  = gmm_standard_mean(X, oldGamma);
  oldPRE = gmm11R_precision(X, oldGamma, oldMU, lambda);
  parS   = 0;
  for (int k=0;k<K;k++){
    for (int i=0;i<P;i++){
      for (int j=i;j<P;j++){
        valijk = std::abs(oldPRE(i,j,k));
        if (valijk >= thres0){
          parS += 1;
        }
      }
    }
  }
  oldPar = par0 + parS;
  for (int k=0;k<K;k++){
    oldSIG.slice(k) = arma::pinv(oldPRE.slice(k));
  }
  if (usediag==true){
    for (int k=0; k<K; k++){
      oldSIG.slice(k) = arma::diagmat(oldSIG.slice(k));
    }
  }
  oldLKD = gmm_loglkd(X, oldPI, oldMU, oldSIG);

  
  // main iteration
  for (int it=0;it<mit;it++){
    // E-step. update Gamma
    newGamma = gmm_standard_gamma(X, oldMU, oldSIG, oldPI);
    // M-step. update parameters + log-likelihood
    newPI  = gmm_standard_pi(newGamma);
    newMU  = gmm_standard_mean(X, newGamma);
    newPRE = gmm11R_precision(X, newGamma, newMU, lambda);
    parS   = 0;
    for (int k=0;k<K;k++){
      for (int i=0;i<P;i++){
        for (int j=i;j<P;j++){
          valijk = std::abs(newPRE(i,j,k));
          if (valijk >= thres0){
            parS += 1;
          }
        }
      }
    }
    newPar = par0 + parS;
    for (int k=0;k<K;k++){
      newSIG.slice(k) = arma::pinv(newPRE.slice(k));
    }
    if (usediag==true){
      for (int k=0; k<K; k++){
        newSIG.slice(k) = arma::diagmat(newSIG.slice(k));
      }
    }
    newLKD = gmm_loglkd(X, newPI, newMU, newSIG);
    // break part
    if ((it>=1)&&(newLKD <= oldLKD)){
      break;
    } else {
      oldGamma = newGamma;
      oldPI    = newPI;
      oldMU    = newMU;
      oldSIG   = newSIG;
      oldLKD   = newLKD;
      oldPRE   = newPRE;
      oldPar   = newPar;
    }
  }
  
  // USE ARMADILLO ROUTINES
  arma::gmm_full model;
  model.reset(P, K);
  model.set_means(arma::trans(oldMU));
  model.set_fcovs(oldSIG);
  model.set_hefts(arma::trans(oldPI));
  
  // RETURN
  Rcpp::List output;
  output["means"]   = oldMU;
  output["covs"]    = oldSIG;
  output["weight"]  = oldPI;
  output["loglkd"]  = model.sum_log_p(arma::trans(X));
  output["cluster"] = arma::trans(model.assign(arma::trans(X), prob_dist));
  return(output);
}

// (03) gmm_16Gfix =============================================================
arma::mat gmm_16Gfix_mean(arma::mat X, arma::mat Gamma, arma::vec Weight){
  // parameters
  int N = Gamma.n_rows;
  int K = Gamma.n_cols;
  int P = X.n_cols;
  
  // iterate
  arma::mat output(K,P,fill::zeros);
  arma::rowvec tmpvec(P,fill::zeros);
  double tmpval=0.0;
  for (int k=0;k<K;k++){
    tmpvec.fill(0.0);
    tmpval = 0.0;
    for (int n=0;n<N;n++){
      tmpvec += Gamma(n,k)*X.row(n)*Weight(n);
      tmpval += Gamma(n,k)*Weight(n);
    }
    output.row(k) = tmpvec/tmpval;
  }
  return(output);
}
arma::cube gmm_16Gfix_cov(arma::mat X, arma::mat Gamma, arma::mat Mu, arma::vec Weight, bool usediag=false){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  int K = Mu.n_rows;
  
  // iterate
  arma::cube output(P,P,K,fill::zeros);
  arma::mat  mslice(P,P,fill::zeros);
  arma::rowvec xdiff(P,fill::zeros);
  arma::mat  tmpmat(P,P,fill::zeros);
  double Nk = 0.0;
  for (int k=0;k<K;k++){
    // denominator
    Nk = arma::accu(Gamma.col(k));
    // numerator
    mslice.fill(0.0);
    for (int n=0;n<N;n++){
      xdiff   = X.row(n) - Mu.row(k);
      mslice += Gamma(n,k)*(xdiff.t()*xdiff)*Weight(n);
    }
    tmpmat = mslice/Nk;
    if (usediag==true){
      output.slice(k) = arma::diagmat(tmpmat);
    } else {
      output.slice(k) = tmpmat;
    }
  }
  return(output);
}
double gmm_16Gfix_loglkd(arma::mat X, arma::vec kPi, arma::mat kMu, arma::cube kSig, arma::vec weight){
  // Parameters
  int n = X.n_rows;
  int p = X.n_cols;
  int k = kMu.n_rows;
  
  // Compute Elementary
  arma::mat piNN(n,k,fill::zeros);
  for (int i=0; i<n; i++){
    for (int j=0; j<k; j++){
      piNN(i,j) = kPi(j)*eval_gaussian_single(X.row(i), kMu.row(j), kSig.slice(j)/weight(i), false);
    }
  }
  
  // Compute
  arma::vec before_log(n,fill::zeros);
  for (int i=0; i<n; i++){
    before_log(i) = std::log(arma::accu(piNN.row(i)));
  }
  
  // Return
  return(arma::accu(before_log));
}
arma::uvec gmm_16Gfix_label(arma::mat& X, arma::mat parMU, arma::cube parSIG, arma::vec parPI, arma::vec weight){
  // parameters
  int N = X.n_rows;
  int K = parSIG.n_slices;
  // compute gamma
  arma::mat piNN(N,K,fill::zeros);
  for (int n=0; n<N; n++){
    for (int k=0; k<K; k++){
      piNN(n,k) = parPI(k)*eval_gaussian_single(X.row(n), parMU.row(k), parSIG.slice(k)/weight(n), false);
    }
  }
  // for each component, find the maximal
  arma::uvec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::index_max(piNN.row(n));
  }
  return(output);
}

// [[Rcpp::export]]
Rcpp::List gmm_16Gfix(arma::mat& X, int k, arma::vec weight, int maxiter, bool usediag){
  // PARAMETERS
  int n = X.n_rows; 
  int p = X.n_cols;

  // PREPARATION : DEFINITION
  arma::mat oldGamma(n,k,fill::zeros);  // (NxK) class membership
  arma::mat newGamma(n,k,fill::zeros);
  arma::mat oldMU(k,p,fill::zeros);     // (KxP) row class means
  arma::mat newMU(k,p,fill::zeros);
  arma::cube oldSIG(p,p,k,fill::zeros); // (PxPxK) covariance slices
  arma::cube newSIG(p,p,k,fill::zeros);
  arma::vec oldPI(k,fill::zeros);       // (K) proportion
  arma::vec newPI(k,fill::zeros);
  double oldLKD = 0.0;
  double newLKD = 0.0;
  
  // PREPARATION : INITIALIZATION
  arma::urowvec initlabel = label_gmm(X, k, 5); // small number of simple task
  for (int i=0; i<n; i++){
    oldGamma(i,initlabel(i)) = 1.0;
  }
  oldPI  = gmm_standard_pi(oldGamma);
  oldMU  = gmm_16Gfix_mean(X, oldGamma, weight);
  oldSIG = gmm_16Gfix_cov(X, oldGamma, oldMU, weight, usediag);
  oldLKD = gmm_16Gfix_loglkd(X, oldPI, oldMU, oldSIG, weight);
  
  // ITERATION
  for (int it=0; it<maxiter; it++){
    // {E} update gamma
    newGamma = gmm_standard_gamma(X, oldMU, oldSIG, oldPI);
    // {M} update parameters + compute loglkd
    newPI  = gmm_standard_pi(newGamma);
    newMU  = gmm_16Gfix_mean(X, newGamma, weight);
    newSIG = gmm_16Gfix_cov(X, newGamma, newMU, weight, usediag);
    newLKD = gmm_16Gfix_loglkd(X, newPI, newMU, newSIG, weight);
    // Breaking Part
    if ((it>0)&&(newLKD <= oldLKD)){
      break;
    } else {
      oldGamma = newGamma;
      oldPI    = newPI;
      oldMU    = newMU;
      oldSIG   = newSIG;
      oldLKD   = newLKD; 
    }
  }
  
  // USE ARMADILLO ROUTINES
  arma::gmm_full model;
  model.reset(p, k);
  model.set_means(arma::trans(oldMU));
  model.set_fcovs(oldSIG);
  model.set_hefts(arma::trans(oldPI));
  
  // RETURN
  Rcpp::List output;
  output["means"]   = oldMU;
  output["covs"]    = oldSIG;
  output["weight"]  = oldPI;
  output["loglkd"]  = model.sum_log_p(arma::trans(X));
  output["cluster"] = gmm_16Gfix_label(X, oldMU, oldSIG, oldPI, weight);
  return(output);
}

// gmm_combine_wsum   : weighted sum -------------------------------------------
// [[Rcpp::export]]
Rcpp::List gmm_combine_wsum(Rcpp::List &gmmlist, arma::vec &weight){
  int K = weight.n_elem;
  
  // access the first one
  Rcpp::List tgtgmm   = gmmlist[0];
  arma::vec  myweight = Rcpp::as<arma::vec>(tgtgmm["weight"])*weight(0);
  arma::mat  mymean   = Rcpp::as<arma::mat>(tgtgmm["mean"]);
  arma::cube myvar    = Rcpp::as<arma::cube>(tgtgmm["variance"]);
  
  // iterate
  arma::vec  tmp_pi;
  arma::mat  tmp_mu;
  arma::cube tmp_var;
  for (int k=1; k<K; k++){
    tgtgmm = gmmlist[k];
    tmp_pi.reset();
    tmp_mu.reset();
    tmp_var.reset();
    
    tmp_pi  = Rcpp::as<arma::vec>(tgtgmm["weight"])*weight(k);
    tmp_mu  = Rcpp::as<arma::mat>(tgtgmm["mean"]);
    tmp_var = Rcpp::as<arma::cube>(tgtgmm["variance"]);
    
    myweight = arma::join_vert(myweight, tmp_pi);
    mymean   = arma::join_vert(mymean,   tmp_mu);
    myvar    = arma::join_slices(myvar,  tmp_var);
  }
  
  Rcpp::List output;
  output["means"]   = mymean;
  output["covs"]    = myvar;
  output["weight"]  = myweight;
  return(output);
}

// gmm_density -----------------------------------------------------------------
// [[Rcpp::export]]
arma::vec gmm_density(arma::mat &coords, arma::vec &weight, arma::mat &mean, arma::cube &variance){
  // parameters
  int N = coords.n_rows;
  int p = coords.n_cols;
  int K = weight.n_elem;
  
  // evaluate per Gaussian + weight
  arma::vec myweight = weight/arma::accu(weight);
  arma::mat eval_class(N,K,fill::zeros);
  for (int k=0; k<K; k++){
    eval_class.col(k) = eval_gaussian(coords, mean.row(k), variance.slice(k), false)*myweight(k);
  }
  
  // finalize
  arma::vec output(N,fill::zeros);
  for (int n=0; n<N; n++){
    output(n) = arma::accu(eval_class.row(n));
  }
  return(output);
}
  
// gmm_pdist_wass2 : compute pairwise distance between GMM components (Wass2) --
// [[Rcpp::export]]
arma::mat gmm_pdist_wass2(arma::mat& mean, arma::cube& variance){
  int N = variance.n_slices;
  int p = variance.n_rows;
  
  arma::cube varsqrt(p,p,N,fill::zeros);
  for (int n=0; n<N; n++){
    varsqrt.slice(n) = arma::sqrtmat_sympd(variance.slice(n));
  }
  
  arma::mat output(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = gauss2dist_wass2(mean.row(i),variance.slice(i),mean.row(j),variance.slice(j),varsqrt.slice(j));
      output(j,i) = output(i,j);
    }
  }
  return(output);
}

// gmm_w2barycenter : W2 barycenter of Gaussian distributions ------------------
// weight : (K)         vector
// mean   : (K x P)     matrix of row-stacked means
// vars   : (P x P x K) cube of stacked covariance matrices

// [[Rcpp::export]]
Rcpp::List gmm_w2barycenter(arma::vec &weight, arma::mat &mean, arma::cube &vars){
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
  Rcpp::List output;
  output["mean"]     = out_mean;
  output["variance"] = oldvar;
  return(output);
}


// (04) gmm_03F : single run of probability ====================================
arma::mat gmm_03F_single(arma::mat Xproj, int kk, int maxiter, bool usediag){
  // parameters
  int nobs = Xproj.n_rows;
  arma::mat Xlowt = arma::trans(Xproj);
  
  // Run a model with case branching
  arma::mat probmat(nobs, kk, fill::zeros);
  if (usediag==true){ // use diagonal
    arma::gmm_diag model;
    bool status = model.learn(Xlowt, kk, maha_dist, random_subset, 5, maxiter, 1e-11, false);
    if (status==false){
      Rcpp::stop("* gmm03F : Fitting GMM with diagonal covariance failed.");
    } else {
      // successful, let's take out the elements and return
      
      arma::mat diagcovs = model.dcovs;
      int p = diagcovs.n_rows;
      int k = diagcovs.n_cols;
      arma::cube myfcovs(p,p,k,fill::zeros);
      for (int i=0; i<k; i++){
        myfcovs.slice(i) = arma::diagmat(diagcovs.col(i));
      }
      probmat = gmm_standard_gamma(Xproj, arma::trans(model.means), myfcovs, arma::trans(model.hefts));
    }
  } else {
    arma::gmm_full model;
    bool status = model.learn(Xlowt, kk, maha_dist, random_subset, 5, maxiter, 1e-11, false);
    if (status==false){
      Rcpp::stop("* gmm03F : Fitting GMM with full covariance failed.");
    } else {
      probmat = gmm_standard_gamma(Xproj, arma::trans(model.means), model.fcovs, arma::trans(model.hefts));
    }
  }
  
  // pairwise similarity
  arma::rowvec tgti(kk,fill::zeros);
  arma::rowvec tgtj(kk,fill::zeros);
  arma::mat simmat(nobs,nobs,fill::zeros);
  for (int i=0; i<(nobs-1); i++){
    tgti = probmat.row(i);
    for (int j=(i+1); j<nobs; j++){
      tgtj = probmat.row(j);
      simmat(i,j) = arma::dot(tgti,tgtj);
      simmat(j,i) = simmat(i,j);
    }
  }
  
  return(simmat);
}

// [[Rcpp::export]]
arma::mat gmm_03F(arma::mat& X, int k, int maxiter, bool usediag, int lowdim, int nruns){
  // parameters
  int n = X.n_rows; 
  int p = X.n_cols;
  arma::mat proj(p,lowdim,arma::fill::zeros);
  arma::mat Xproj(n,lowdim,arma::fill::zeros);
  
  arma::cube simcube(n,n,nruns,fill::zeros);
  arma::mat  simmat(n,n,fill::zeros);
  
  
  // iterations
  for (int it=0; it<nruns; it++){
    // random projection
    proj.randn();
    Xproj = X*proj; 
    
    // compute similarity matrix
    simcube.slice(it) = gmm_03F_single(Xproj, k, maxiter, usediag);
  }
  
  // return
  simmat = arma::mean(simcube, 2);
  return(simmat);
}