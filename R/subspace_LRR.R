#' Low-Rank Representation
#' 
#' Low-Rank Representation (LRR) constructs the connectivity of the data by 
#' solving 
#' \deqn{\textrm{min}_C \|C\|_*\quad\textrm{such that}\quad D=DC}
#' for column-stacked data matrix \eqn{D} and \eqn{\|\cdot \|_*} is the 
#' nuclear norm which is relaxation of the rank condition. If you are interested in 
#' full implementation of the algorithm with sparse outliers and noise, please 
#' contact the maintainer.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param rank sum of dimensions for all \eqn{k} subspaces (default: 2).
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
#' ## run LRR algorithm with k=2, 3, and 4 with rank=4
#' output2 = LRR(data, k=2, rank=4)
#' output3 = LRR(data, k=3, rank=4)
#' output4 = LRR(data, k=4, rank=4)
#' 
#' ## extract label information
#' lab2 = output2$cluster
#' lab3 = output3$cluster
#' lab4 = output4$cluster
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(dat2, pch=19, cex=0.9, col=lab2, main="LRR:K=2")
#' plot(dat2, pch=19, cex=0.9, col=lab3, main="LRR:K=3")
#' plot(dat2, pch=19, cex=0.9, col=lab4, main="LRR:K=4")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{liu_robust_2010}{T4cluster}
#' 
#' @concept subspace
#' @export
LRR <- function(data, k=2, rank=2){
  ## PREPARE : EXPLICIT INPUTS
  X = prec_input_matrix(data)
  myk = max(1, round(k))
  myr = max(1, round(rank))
  
  ## COMPUTE EVERYTHING IN C++
  cpprun = cpp_LRR(X, myk, myr)
  
  ## WRAP
  output  = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  # output$eigval  = as.vector(cpprun$values)
  # output$embeds  = cpprun$embeds
  output$algorithm = "LRR"
  return(structure(output, class="T4cluster"))
}
