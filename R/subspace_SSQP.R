#' Subspace Segmentation via Quadratic Programming
#' 
#' Subspace Segmentation via Quadratic Programming (SSQP) solves the following problem
#' \deqn{\textrm{min}_Z \|X-XZ\|_F^2 + \lambda \|Z^\top Z\|_1 \textrm{ such that }diag(Z)=0,~Z\leq 0}
#' where \eqn{X\in\mathbf{R}^{p\times n}} is a column-stacked data matrix. The computed \eqn{Z^*} is 
#' used as an affinity matrix for spectral clustering.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param lambda regularization parameter (default: 1e-5).
#' @param ... extra parameters for the gradient descent algorithm including \describe{
#' \item{maxiter}{maximum number of iterations (default: 100).}
#' \item{abstol}{tolerance level to stop (default: 1e-7).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## generate a toy example
#' set.seed(10)
#' tester = genLP(n=100, nl=2, np=1, iso.var=0.1)
#' data   = tester$data
#' label  = tester$class
#' 
#' ## do PCA for data reduction
#' proj = base::eigen(stats::cov(data))$vectors[,1:2]
#' dat2 = data%*%proj
#' 
#' ## run SSQP for k=3 with different lambda values
#' out1 = SSQP(data, k=3, lambda=1e-2)
#' out2 = SSQP(data, k=3, lambda=1)
#' out3 = SSQP(data, k=3, lambda=1e+2)
#' 
#' ## extract label information
#' lab1 = out1$cluster
#' lab2 = out2$cluster
#' lab3 = out3$cluster
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(dat2, pch=19, cex=0.9, col=lab1, main="SSQP:lambda=1e-2")
#' plot(dat2, pch=19, cex=0.9, col=lab2, main="SSQP:lambda=1")
#' plot(dat2, pch=19, cex=0.9, col=lab3, main="SSQP:lambda=1e+2")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{wang_efficient_2011}{T4cluster}
#' 
#' @concept subspace
#' @export
SSQP <- function(data, k=2, lambda=1e-5, ...){
  ## PRELIMINARY
  X = prec_input_matrix(data)
  K = max(1, round(k))
  mylbd  = max(10*.Machine$double.eps, as.double(lambda))
  
  params = list(...)
  pnames = names(params)
  myiter = ifelse(("maxiter"%in%pnames), params$maxiter, 100)
  myiter = max(5, round(myiter))
  mytol  = ifelse(("abstol"%in%pnames), params$abstol, 1e-7)
  mytol  = max(as.double(mytol), sqrt(.Machine$double.eps))
  
  ## RUN C++
  cpprun = cpp_SSQP(X, K, mylbd, myiter, mytol)
  
  ## WRAP
  output = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  output$algorithm = "SSQP"
  return(structure(output, class="T4cluster"))
}