#' Spectral Clustering by Gu and Wang (2009)
#' 
#' The algorithm defines a set of data-driven 
#' bandwidth parameters where \eqn{\sigma_i} is the average distance from a point \eqn{x_i} to its \code{nnbd}-th 
#' nearest neighbor. Then the affinity matrix is defined as
#' \deqn{A_{ij} = \exp(-d(x_i, d_j)^2 / \sigma_i \sigma_j)} and the standard 
#' spectral clustering of Ng, Jordan, and Weiss (\code{\link{scNJW}}) is applied.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations or S3 \code{dist} object of \eqn{n} observations.
#' @param k the number of clusters (default: 2).
#' @param nnbd neighborhood size to define data-driven bandwidth parameter (default: 7).
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
#' cl2 = sc09G(X, k=2)$cluster
#' cl3 = sc09G(X, k=3)$cluster
#' cl4 = sc09G(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="sc09G: k=2")
#' plot(X2d, col=cl3, pch=19, main="sc09G: k=3")
#' plot(X2d, col=cl4, pch=19, main="sc09G: k=4")
#' par(opar)
#' 
#' @references 
#' \insertRef{gu_improved_2009}{T4cluster}
#' 
#' @concept algorithm
#' @export
sc09G <- function(data, k=2, nnbd=7, ...){
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_dist(data)
  myk     = max(1, round(k))
  mynnbd  = max(5, round(nnbd))
  
  ## PREPARE : IMPLICIT ONES
  params  = list(...)
  pnames  = names(params)
  
  myiter     = ifelse(("maxiter"%in%pnames), max(10, round(params$maxiter)), 10)
  myclust    = ifelse(("algclust"%in%pnames), match.arg(tolower(params$algclust),c("kmeans","gmm")), "kmeans") 
  kmeansflag = ifelse(all(myclust=="kmeans"), TRUE, FALSE)
  
  ## RUN
  cpprun = cpp_sc09G(mydata, myk, mynnbd, kmeansflag, myiter)
  
  ## WRAP
  output  = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  output$eigval  = as.vector(cpprun$values)
  output$embeds  = cpprun$embeds
  output$algorithm = "sc09G"
  return(structure(output, class="T4cluster"))
}