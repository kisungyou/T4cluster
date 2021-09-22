#' Generalized Bayesian Clustering with PHATE Geometry
#' 
#' 
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param ... extra parameters including \describe{
#' \item{nnbd}{size of nearest neighborhood (default: 5).}
#' \item{alpha}{decay parameter for Gaussian kernel exponent (default: 2).}
#' \item{gamma_a}{hyperparameter of the Gamma prior (default: 0).}
#' \item{gamma_b}{hyperparameter of the Gamma prior (default: 0).}
#' \item{burn_in}{number of MCMC samples to be thrown away for burn-in (default: 50).}
#' \item{nsample}{number of MCMC samples to be kept after burn-in (default: 50).}
#' \item{random.start}{a logical to use random initialization or not (default: \code{TRUE}).}
#' \item{print.progress}{a logical to show how much an algorithm is proceeded (default: \code{FALSE}.)}
#' }
#' 
#' 
#' @concept algorithm
#' @export
gbphate <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUTS
  mydata = prec_input_matrix(data)
  myk    = max(1, round(k))
  myn    = base::nrow(mydata)
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  if ("nnbd"%in%pnames){
    my_nnbd = max(1, round(pars$nnbd))
  } else {
    my_nnbd = 5
  }
  if ("alpha"%in%pnames){
    my_alpha = max(1, as.double(pars$alpha))
  } else {
    my_alpha = 2
  }
  if ("gamma_a"%in%pnames){
    my_gamma_a = max(0, as.double(pars$gamma_a))
  } else {
    my_gamma_a = 0.0
  }
  if ("gamma_b"%in%pnames){
    my_gamma_b = max(0, as.double(pars$gamma_b))
  } else {
    my_gamma_b = 0.0
  }
  if ("burn_in"%in%pnames){
    my_burn_in = max(5, round(pars$burn_in))
  } else {
    my_burn_in = 50
  }
  if ("nsample"%in%pnames){
    my_nsample = max(5, round(pars$nsample))
  } else {
    my_nsample = 50
  }
  if ("random.start"%in%pnames){
    my_randstart = as.logical(pars$random.start)
  } else {
    my_randstart = TRUE
  }
  
  
  print.progress = TRUE
  
  ## RUN PHATE FOR OPTIMAL ROW-STOCHASTIC MATRIX RETRIEVAL
  if ((is.matrix(mydata))&&all(rowSums(mydata)==1)){
    pseudoX <- base::sqrt(mydata)
  } else {
    distX <- stats::dist(mydata)
    fun_phate <- utils::getFromNamespace("hidden_PHATE","maotai")
    run_phate <- fun_phate(distX, nbdk=my_nnbd, alpha=my_alpha)
    pseudoX   <- base::sqrt(run_phate$P)
  }
  if (print.progress){
    print("* gbphate : Stage 1 for PHATE geometry recovery is complete.")
  }
  
  ## PREPARE FOR RUNNING WITH GBCLUST
  if (my_randstart){
    map_kmeans <- base::sample(1:myk, myn, replace = TRUE)
  } else {
    map_kmeans <- kmeans_fast(pseudoX, myk)
  }
  map_freq   <- as.vector(base::table(map_kmeans))
  
  ## FIT AND RUN THE GIBBS CLUSTERING
  gbrun <- spkmeans_gibbs(my_burn_in, my_nsample, pseudoX,
                          my_gamma_a, my_gamma_b,
                          map_kmeans, map_freq, print.progress) # returns lambda estimates
  
  ## REARRANGE FOR RETURN
  output = list()
  output$pseudoX  = pseudoX
  output$clusters = round(gbrun$G)
  output$lambdas  = as.vector(gbrun$lambda)
  output$losses   = as.vector(gbrun$loss) 
  output$algorithm ="gbphate"
  
  if (print.progress){
    print("* gbphate : Stage 3 is done. Algorithm terminated.")
  }
  return(structure(output, class="T4cluster"))
}



