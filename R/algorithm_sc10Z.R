#' Spectral Clustering by Zhang et al. (2010)
#' 
#' The algorithm defines a set of data-driven 
#' bandwidth parameters \eqn{p_{ij}} by constructing a similarity matrix. 
#' Then the affinity matrix is defined as \deqn{A_{ij} = \exp(-d(x_i, d_j)^2 / 2 p_{ij})} 
#' and the standard spectral clustering of Ng, Jordan, and Weiss (\code{\link{scNJW}}) is applied.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations or S3 \code{dist} object of \eqn{n} observations.
#' @param k the number of clusters (default: 2).
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
#' cl2 = sc10Z(X, k=2)$cluster
#' cl3 = sc10Z(X, k=3)$cluster
#' cl4 = sc10Z(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="sc10Z: k=2")
#' plot(X2d, col=cl3, pch=19, main="sc10Z: k=3")
#' plot(X2d, col=cl4, pch=19, main="sc10Z: k=4")
#' par(opar)
#' 
#' @references 
#' \insertRef{zhang_spectral_2010}{T4cluster}
#' 
#' @concept algorithm
#' @export
sc10Z <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_dist(data)
  myk     = max(1, round(k))

  ## PREPARE : IMPLICIT ONES
  params  = list(...)
  pnames  = names(params)
  
  myiter     = ifelse(("maxiter"%in%pnames), max(10, round(params$maxiter)), 10)
  myclust    = ifelse(("algclust"%in%pnames), match.arg(tolower(params$algclust),c("kmeans","gmm")), "kmeans") 
  kmeansflag = ifelse(all(myclust=="kmeans"), TRUE, FALSE)
  
  ## RUN
  cpprun = cpp_sc10Z(mydata, myk, kmeansflag, myiter)
  
  ## WRAP
  output  = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  output$eigval  = as.vector(cpprun$values)
  output$embeds  = cpprun$embeds
  output$algorithm = "sc10Z"
  return(structure(output, class="T4cluster"))
}