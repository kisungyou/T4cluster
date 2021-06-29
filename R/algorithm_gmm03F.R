#' Ensemble of Gaussian Mixtures with Random Projection
#' 
#' When the data lies in a high-dimensional Euclidean space, fitting a model-based 
#' clustering algorithm is troublesome. This function implements an algorithm 
#' from the reference, which uses an aggregate information from an ensemble of 
#' Gaussian mixtures in combination with random projection.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param ... extra parameters including \describe{
#' \item{nruns}{the number of projections (default: 20).}
#' \item{lowdim}{target dimension for random projection (default: 5).}
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{usediag}{a logical; covariances are diagonal if \code{TRUE}, or full covariances are returned for \code{FALSE} (default: \code{FALSE}).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).}
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @references 
#' \insertRef{10.5555/3041838.3041862}{T4cluster}
#' 
#' @examples
#' \donttest{
#' # -------------------------------------------------------------
#' #            clustering with 'iris' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(iris)
#' X   = as.matrix(iris[,1:4])
#' lab = as.integer(as.factor(iris[,5]))
#' 
#' ## EMBEDDING WITH PCA
#' X2d = Rdimtools::do.pca(X, ndim=2)$Y  
#' 
#' ## CLUSTERING WITH DIFFERENT K VALUES
#' cl2 = gmm03F(X, k=2)$cluster
#' cl3 = gmm03F(X, k=3)$cluster
#' cl4 = gmm03F(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="gmm03F: k=2")
#' plot(X2d, col=cl3, pch=19, main="gmm03F: k=3")
#' plot(X2d, col=cl4, pch=19, main="gmm03F: k=4")
#' par(opar)
#' }
#' 
#' @concept algorithm
#' @export
gmm03F <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_matrix(data)
  myk     = max(1, round(k))
  myn     = base::nrow(data)
  myp     = base::ncol(data)
  if (myp < 2){
    stop("* gmm03F : for univariate data, use other functions.")
  }
  
  ## PREPARE : IMPLICIT
  params  = list(...)
  pnames  = names(params)
  
  myiter = round(max(5, ifelse(("maxiter"%in%pnames), params$maxiter, 10)))
  mydiag = as.logical(ifelse(("usediag"%in%pnames), params$usediag, FALSE))
  
  if ("lowdim"%in%pnames){
    mylowdim = max(2,round(params$lowdim))
  } else {
    mylowdim = 5
  }
  if (mylowdim > myp){
    mylowdim = 2
  }
  
  if ("nruns"%in%pnames){
    mynruns = max(5, round(params$nruns))
  } else {
    mynruns = 20
  }
  
  ## MAIN COMPUTATION FOR SIMILARITY
  mat.sim = gmm_03F(mydata, myk, myiter, mydiag, mylowdim, mynruns)
  mat.P   = 1-mat.sim; diag(mat.P) = 0; mat.P[(mat.P<0)] = 0

  ## HCLUST+COMPLETE LINKAGE
  fimport = utils::getFromNamespace("hidden_hclust", "maotai")
  hcout   = fimport(stats::as.dist(mat.P), "complete", NULL)
  
  ## WRAP AND RETURN
  output = list()
  output$cluster  = as.vector(stats::cutree(hcout, k=myk))
  output$algorithm = "gmm03F"
  return(structure(output, class="T4cluster"))  
}
