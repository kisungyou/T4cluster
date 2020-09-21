#' Generate Data from Fitted GMM Models
#' 
#' Given a fitted GMM model in \eqn{\mathbb{R}^p}, it generates \eqn{n} samples 
#' according to the fitted model. 
#' 
#' @param n the number of observations to be drawn.
#' @param gmmobj an output of any GMM routines in our package of \code{T4cluster} class.
#' 
#' @return an \eqn{(n\times p)} generated data matrix.
#' 
#' @examples
#' # -------------------------------------------------------------
#' #            clustering with 'iris' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(iris)
#' X   = as.matrix(iris[,1:4])
#' 
#' ## FIT THE MODEL WITH K=3
#' cl3 = gmm(X, k=4)
#' 
#' ## GENERATE 150 SAMPLES FROM THE MODEL
#' Z = gmmsample(150, cl3)
#' 
#' ## PCA FOR ORIGINAL AND GENERATED DATA
#' X2trf = Rdimtools::do.pca(X, ndim=2)
#' X2ori = X2trf$Y
#' X2gen = (Z-matrix(rep(colMeans(Z),150),ncol=4,byrow=TRUE))%*%X2trf$projection
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(X2ori, pch=19, cex=0.8, main="original data")
#' plot(X2gen, pch=19, cex=0.8, main="generated data")
#' par(opar)
#' 
#' @concept utility
#' @export
gmmsample <- function(n, gmmobj){
  ## PREPARE : EXPLICIT INPUT
  myn = max(1, round(n))
  if (!inherits(gmmobj,"T4cluster")){
    stop("* gmmsample : input 'gmmobj' should be an object from any gmm functions in our package.")
  }
  if ((!("algorithm"%in%names(gmmobj)))||(!(all(gmm_name_segment(gmmobj$algorithm)=="gmm")))){
    stop("* gmmsample : the given 'gmmobj' is not an output of gmm functions in our package.")
  }
  
  ## RUN
  output = gmm_sample(myn, gmmobj$weight, gmmobj$mean, gmmobj$variance)
  return(output)
}