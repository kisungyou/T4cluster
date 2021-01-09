#' Spectral Clustering by Yang et al. (2011)
#' 
#' As a data-driven method, the algorithm recovers geodesic distance from a k-nearest 
#' neighbor graph scaled by an (exponential) parameter \eqn{\rho} and applies 
#' random-walk spectral clustering. Authors referred their method as 
#' density sensitive similarity function.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations or S3 \code{dist} object of \eqn{n} observations.
#' @param k the number of clusters (default: 2).
#' @param nnbd neighborhood size to define data-driven bandwidth parameter (default: 7).
#' @param rho exponent scaling parameter (default: 2).
#' @param ... extra parameters including \describe{
#' \item{algclust}{method to perform clustering on embedded data; either \code{"kmeans"} (default) or \code{"GMM"}.}
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{eigval}{eigenvalues of the graph laplacian's spectral decomposition.}
#' \item{embeds}{an \eqn{(n\times k)} low-dimensional embedding.}
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @examples
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
#' cl2 = sc11Y(X, k=2)$cluster
#' cl3 = sc11Y(X, k=3)$cluster
#' cl4 = sc11Y(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="sc11Y: k=2")
#' plot(X2d, col=cl3, pch=19, main="sc11Y: k=3")
#' plot(X2d, col=cl4, pch=19, main="sc11Y: k=4")
#' par(opar)
#' 
#' @references 
#' \insertRef{yang_spectral_2011}{T4cluster}
#' 
#' @concept algorithm
#' @export
sc11Y <- function(data, k=2, nnbd=7, rho=2, ...){
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_dist(data)
  myk     = max(1, round(k))
  mynbd   = max(7, round(nnbd))
  myrho   = max(1, as.double(rho))
  
  ## PREPARE : IMPLICIT ONES
  params  = list(...)
  pnames  = names(params)
  
  myiter     = ifelse(("maxiter"%in%pnames), max(10, round(params$maxiter)), 10)
  myclust    = ifelse(("algclust"%in%pnames), match.arg(tolower(params$algclust),c("kmeans","gmm")), "kmeans") 
  kmeansflag = ifelse(all(myclust=="kmeans"), TRUE, FALSE)
  
  ## MAIN COMPUTATION
  # 1. nearest neighbor search
  nbdfunc  = utils::getFromNamespace("hidden_knn","maotai")
  ksearch  = nbdfunc(mydata, nnbd=(mynbd+1))
  # 2. run 'cpp_sc11Y'
  cpprun = cpp_sc11Y(ksearch$nn.idx[,2:(mynbd+1)], ksearch$nn.dists[,2:(mynbd+1)], myk, kmeansflag, myiter, myrho)
  
  # Wrap and Return
  output  = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  output$eigval  = as.vector(cpprun$values)
  output$embeds  = cpprun$embeds
  output$algorithm = "sc11Y"
  return(structure(output, class="T4cluster"))  
}