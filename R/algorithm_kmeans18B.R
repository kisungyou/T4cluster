#' K-Means Clustering with Lightweight Coreset
#' 
#' Apply \eqn{k}-means clustering algorithm on top of the lightweight coreset 
#' as proposed in the paper. 
#' The smaller the set is, the faster the execution becomes with potentially larger quantization errors.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param m the size of coreset (default: \eqn{n/2}).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{nstart}{the number of random initializations (default: 5).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{mean}{a \eqn{(k\times p)} matrix where each row is a class mean.}
#' \item{wcss}{within-cluster sum of squares (WCSS).}
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
#' ## CLUSTERING WITH DIFFERENT CORESET SIZES WITH K=3
#' core1 = kmeans18B(X, k=3, m=25)$cluster
#' core2 = kmeans18B(X, k=3, m=50)$cluster
#' core3 = kmeans18B(X, k=3, m=100)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=core1, pch=19, main="kmeans18B: m=25")
#' plot(X2d, col=core2, pch=19, main="kmeans18B: m=50")
#' plot(X2d, col=core3, pch=19, main="kmeans18B: m=100")
#' par(opar)
#' 
#' @references 
#' \insertRef{bachem_scalable_2018}{T4cluster}
#' 
#' @concept algorithm
#' @export
kmeans18B <- function(data, k=2, m=round(nrow(data)/2), ...){
  ## PREPARE : EXPLICIT INPUTS
  mydata = prec_input_matrix(data)
  myn    = base::nrow(mydata)
  myk    = max(1, round(k))
  mym    = max(2*myk, round(m))
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  myiter   = ifelse(("maxiter" %in% pnames), max(10, round(pars$maxiter)), 10)
  mynstart = ifelse(("nstart"%in%pnames), max(1, round(pars$nstart)), 5)
  
  ## MULTIPLE RUNS
  multiple_runs = list()
  multiple_wcss = rep(0, mynstart)
  for (i in 1:mynstart){
    multiple_runs[[i]] = coreset_18B(mydata, myk, mym, myiter)
    multiple_wcss[i]   = multiple_runs[[i]]$wcss
  }
  
  ## SELECT THE BEST RUN
  cpprun = multiple_runs[[which.min(multiple_wcss)]]
  
  ## WRAP AND RETURN
  cpprun$cluster = as.vector(cpprun$cluster)+1
  cpprun$algorithm ="kmeans18B"
  return(structure(cpprun, class="T4cluster"))
}