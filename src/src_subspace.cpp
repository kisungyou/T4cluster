#include "RcppArmadillo.h"
#include "utilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/* SRC FUNCTIONS FOR SUBSPACE CLUSTERING
 * subspace_normalNJW : spectral clustering without zero-ing the diagonal
 * subspace_normalSM
 * 
 * (01) fast_loss_prj : compute loss function for MSM
 * (02) cpp_LRR       : LRR  simplest
 * (03) cpp_LRSC      : LRSC simplest
 * (04) cpp_EKSS_0    : zero iterations
 *      cpp_EKSS_T
 *      cpp_EKSS_affinity
 * (05) cpp_LSR
 * (06) cpp_SSQP
 */


Rcpp::List subspace_normalNJW(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian 
  arma::mat A = W; 
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
Rcpp::List subspace_normalSM(arma::mat W, int K, bool usekmeans, int maxiter){
  // build laplacian
  arma::mat A = W; 
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
  // double Wmax = 0.0;
  
  // compute affinity
  arma::mat Vr = V.head_cols(par_r);
  arma::mat C  = Vr*Vr.t();
  arma::mat W  = arma::abs(C);
  // W.diag().fill(0.0);
  
  // Wmax = W.max();
  // if (Wmax > 1.0){
  //   W /= Wmax;
  // }
  
  
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
  Rcpp::List output = subspace_normalNJW(W, par_k, my_kmeans, my_maxiter);
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
  // W.diag().fill(0.0);
  // double Wmax = W.max();
  // if (Wmax > 1.0){
  //   W /= Wmax;
  // }
  
  bool my_kmeans = true;
  int my_maxiter = 100;
  Rcpp::List output = subspace_normalNJW(W, K, my_kmeans, my_maxiter);
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
    for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
        if (label_now(i)==label_now(j)){
          A(i,j) = A(i,j) + (1.0/BB);
        }
      }
    }
    // for (int i=0; i<(N-1); i++){
    //   for (int j=(i+1); j<N; j++){
    //     if (label_now(i)==label_now(j)){
    //       A(i,j) = A(i,j) + (1.0/BB);
    //       A(j,i) = A(j,i) + (1.0/BB);
    //     }
    //   }
    // }
  }
  return(A);
}

// (05) cpp_LSR ================================================================
// input is directly given in (n x p) form
// [[Rcpp::export]]
Rcpp::List cpp_LSR(arma::mat& data, int K, double lambda, bool zerodiag){
  // preliminary
  int N = data.n_rows;
  arma::mat D = arma::inv(data*data.t() + lambda*arma::eye(N, N)); // (X'X + lbd*I)^{-1}
  arma::mat Z(N,N,fill::zeros);
  
  // computing affinity matrix
  if (zerodiag){
    arma::vec vecDinv = 1.0/arma::diagvec(D);
    Z = -D*arma::diagmat(vecDinv);
    Z.diag().fill(0.0);
  } else {
    Z = D*(data*data.t());
  }
  arma::mat affinity = (arma::abs(Z)+arma::abs(Z.t()))/2.0;
  // if (affinity.max() > 1.0){
  //   affinity /= affinity.max();
  // }
  
  // run NJW spectral clustering & return
  bool my_kmeans = true;
  int my_maxiter = 100;
  Rcpp::List output = subspace_normalSM(affinity, K, my_kmeans, my_maxiter);
  return(output);
}

// (06) cpp_SSQP ===============================================================
// [[Rcpp::export]]
Rcpp::List cpp_SSQP(arma::mat& data, int K, double lambda, int maxiter, double tolerance){
  // preliminary
  arma::mat X = data.t();
  int N = X.n_cols;
  
  // initialization
  arma::mat XtX  = X.t()*X;
  arma::mat Zold = arma::inv(XtX + lambda*arma::eye(N,N))*(X.t()*X);
  Zold.diag().fill(0.0);
  Zold.elem(arma::find(Zold<0.0)).zeros();

  arma::mat Ztmp(N,N,fill::zeros);
  arma::mat Znew(N,N,fill::zeros);
  arma::mat E(N,N,fill::ones);
  
  arma::mat Dtmp(N,N,fill::zeros);
  arma::mat Dfin(N,N,fill::zeros);

  arma::mat grad_old = 2.0*XtX*Zold - 2.0*XtX + 2.0*lambda*Zold*E;
  arma::mat grad_new(N,N,fill::zeros);
  double fval_old = std::pow(arma::norm(X-X*Zold,"fro"), 2.0) + lambda*arma::accu(arma::abs(Zold.t()*Zold));
  double fval_new = 0.0;
  
  int nsteps = 50;
  arma::vec update_step(nsteps,fill::zeros);
  arma::vec update_fval(nsteps,fill::zeros);
  update_step(0) = 2.5;
  for (int i=0; i<(nsteps-1); i++){
    update_step(i+1) = update_step(i)*0.5;
  }
  
  // iteration with Spectral Project Gradient method (1999, Birgin)
  double sigma  = 1.0;
  double incthr = 0.0;
  arma::vec s;
  arma::vec y;
  for (int it=0; it<maxiter; it++){
    // compute D
    Dtmp = Zold - sigma*grad_old;
    Dtmp.diag().fill(0.0);
    Dtmp.elem(arma::find(Dtmp<0.0)).zeros();
    Dfin = Dtmp - Zold;
    
    // optimal updating
    // if no further search object is better, stop
    update_fval.fill(arma::datum::inf);
    for (int i=0; i<nsteps; i++){
      Ztmp = Zold + update_step(i)*Dfin;
      update_fval(i) = std::pow(arma::norm(X-X*Ztmp,"fro"), 2.0) + lambda*arma::accu(arma::abs(Ztmp.t()*Ztmp));
      if (update_fval(i) < fval_old){
        // Rcpp::Rcout << "SSQP : iteration " << it << " - only first " << i << " used." << std::endl;
        break;
      }
    }
    
    fval_new = update_fval.min();
    Znew     = Zold + update_step(update_fval.index_min())*Dfin;
    grad_new = 2.0*XtX*Znew - 2.0*XtX + 2.0*lambda*Znew*E;
    
    // auxiliary information
    s = arma::vectorise(Znew-Zold);
    y = arma::vectorise(grad_new-grad_old);
    sigma = arma::dot(y,y)/arma::dot(s,y);
    
    // updating information
    incthr = std::pow(arma::norm(Zold-Znew,"fro"),2.0);
    Zold = Znew;
    grad_old = grad_new;
    fval_old = fval_new;
    if ((incthr < tolerance)){
      // Rcpp::Rcout << "SSQP : iteration " << it << " breaks!" << std::endl;
      break;
    }
    // Rcpp::Rcout << "SSQP : iteration " << it << " complete..." << std::endl;
  }
  
  // run spectral clustering
  arma::mat affinity = (arma::abs(Zold) + arma::abs(Zold.t()))/2.0;
  bool my_kmeans = true;
  int my_maxiter = 100;
  Rcpp::List output = subspace_normalNJW(affinity, K, my_kmeans, my_maxiter);
  return(output);
}