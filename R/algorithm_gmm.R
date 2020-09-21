#' Finite Gaussian Mixture Model
#' 
#' Finite Gaussian Mixture Model (GMM) is a well-known probabilistic clustering algorithm by fitting the following distribution to the data
#' \deqn{f(x; \left\lbrace \mu_k, \Sigma_k \right\rbrace_{k=1}^K) = \sum_{k=1}^K w_k N(x; \mu_k, \Sigma_k)}
#' with parameters \eqn{w_k}'s for cluster weights, \eqn{\mu_k}'s for class means, and \eqn{\Sigma_k}'s for class covariances. 
#' This function is a wrapper for \pkg{Armadillo}'s GMM function, which supports two types of covariance models.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{usediag}{a logical; covariances are diagonal if \code{TRUE}, or full covariances are returned for \code{FALSE} (default: \code{FALSE}).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).}
#' \item{mean}{a \eqn{(k\times p)} matrix where each row is a class mean.}
#' \item{variance}{a \eqn{(p\times p\times k)} array where each slice is a class covariance.}
#' \item{weight}{a length-\eqn{k} vector of class weights that sum to 1.}
#' \item{loglkd}{log-likelihood of the data for the fitted model.}
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
#' cl2 = gmm(X, k=2)$cluster
#' cl3 = gmm(X, k=3)$cluster
#' cl4 = gmm(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="gmm: k=2")
#' plot(X2d, col=cl3, pch=19, main="gmm: k=3")
#' plot(X2d, col=cl4, pch=19, main="gmm: k=4")
#' par(opar)
#' 
#' @concept algorithm
#' @export
gmm <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_matrix(data)
  myk     = max(1, round(k))
  
  ## PREPARE : IMPLICIT
  params  = list(...)
  pnames  = names(params)
  
  myiter = round(max(5, ifelse(("maxiter"%in%pnames), params$maxiter, 10)))
  mydiag = as.logical(ifelse(("usediag"%in%pnames), params$usediag, FALSE))
  
  ## RUN
  cpprun = gmm_armadillo(data, myk, myiter, mydiag)
  
  ## WRAP AND RETURN
  output = list()
  output$cluster  = as.vector(cpprun$cluster+1)
  output$mean     = cpprun$means
  output$variance = cpprun$covs
  output$weight   = as.vector(cpprun$weight)
  output$loglkd   = as.double(cpprun$loglkd)
  output$algorithm = "gmm"
  return(structure(output, class="T4cluster"))  
}