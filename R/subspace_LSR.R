#' Least Squares Regression
#' 
#' For the subspace clustering, traditional method of least squares regression 
#' is used to build coefficient matrix that reconstructs the data point by solving 
#' \deqn{\textrm{min}_Z \|X-XZ\|_F^2 + \lambda \|Z\|_F \textrm{ such that }diag(Z)=0}
#' where \eqn{X\in\mathbf{R}^{p\times n}} is a column-stacked data matrix. 
#' As seen from the equation, we use a denoising version controlled by \eqn{\lambda} and 
#' provide an option to abide by the constraint \eqn{diag(Z)=0} by \code{zerodiag} parameter. 
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param lambda regularization parameter (default: 1e-5).
#' @param zerodiag a logical; \code{TRUE} (default) to use the problem formulation with zero 
#' diagonal entries or \code{FALSE} otherwise.
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @references 
#' \insertRef{hutchison_robust_2012}{T4cluster}
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
#' ## run LSR for k=3 with different lambda values
#' out1 = LSR(data, k=3, lambda=1e-2)
#' out2 = LSR(data, k=3, lambda=1)
#' out3 = LSR(data, k=3, lambda=1e+2)
#' 
#' ## extract label information
#' lab1 = out1$cluster
#' lab2 = out2$cluster
#' lab3 = out3$cluster
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(dat2, pch=19, cex=0.9, col=lab1, main="LSR:lambda=1e-2")
#' plot(dat2, pch=19, cex=0.9, col=lab2, main="LSR:lambda=1")
#' plot(dat2, pch=19, cex=0.9, col=lab3, main="LSR:lambda=1e+2")
#' par(opar)
#' }
#' 
#' @concept subspace
#' @export
LSR <- function(data, k=2, lambda=1e-5, zerodiag=TRUE){
  ## PREPARE : EXPLICIT INPUTS
  X = prec_input_matrix(data)
  K = max(1, round(k))
  mylbd   = max(10*.Machine$double.eps, as.double(lambda))
  myzdiag = as.logical(zerodiag)
  
  ## RUN C++
  cpprun = cpp_LSR(X, K, mylbd, myzdiag)

  ## WRAP
  output = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  output$algorithm = "LSR"
  return(structure(output, class="T4cluster"))
}