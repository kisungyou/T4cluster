#' DP-Means Clustering
#' 
#' DP-means is a non-parametric clustering method motivated by DP mixture model in that 
#' the number of clusters is determined by a parameter \eqn{\lambda}. The larger 
#' the \eqn{\lambda} value is, the smaller the number of clusters is attained. 
#' In addition to the original paper, we added an option to randomly permute 
#' an order of updating for each observation's membership as a common 
#' heuristic in the literature of cluster analysis. 
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param lambda a threshold to define a new cluster (default: 0.1).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{eps}{the stopping criterion for iterations (default: 1e-5).}
#' \item{permute}{a logical; \code{TRUE} if random order for update is used, \code{FALSE} otherwise (default).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
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
#' ## CLUSTERING WITH DIFFERENT LAMBDA VALUES
#' dpm1 = dpmeans(X, lambda=1)$cluster
#' dpm2 = dpmeans(X, lambda=5)$cluster
#' dpm3 = dpmeans(X, lambda=25)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=dpm1, pch=19, main="dpmeans: lambda=1")
#' plot(X2d, col=dpm2, pch=19, main="dpmeans: lambda=5")
#' plot(X2d, col=dpm3, pch=19, main="dpmeans: lambda=25")
#' par(opar)
#' 
#' @references 
#' \insertRef{kulis_revisiting_2012}{T4cluster}
#' 
#' @concept algorithm
#' @export
dpmeans <- function(data, lambda=0.1, ...){
  ## PREPARE : EXPLICIT INPUTS
  mydata = prec_input_matrix(data)
  mylbd  = max(sqrt(.Machine$double.eps), as.double(lambda))
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  myiter   = max(5, ifelse(("maxiter" %in% pnames), round(pars$maxiter), 10))
  myeps    = max(sqrt(.Machine$double.eps), ifelse(("eps"%in%pnames), pars$eps, 1e-5))
  myperm   = as.logical(ifelse(("permute"%in%pnames), pars$permute, FALSE))
  
  ## RUN DPMEANS ALGORITHM
  dpcluster = maotai::dpmeans(mydata, lambda=mylbd, maxiter=myiter, abstol=myeps, permute.order=myperm)$cluster
  
  ## WRAP AND RETURN
  output = list()
  output$cluster = dpcluster
  output$algorithm ="dpmeans"
  return(structure(output, class="T4cluster"))
}