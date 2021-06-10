#' Sparse Subspace Clustering
#' 
#' Sparse Subspace Clustering (SSC) assumes that the data points lie in 
#' a union of low-dimensional subspaces. The algorithm constructs local 
#' connectivity and uses the information for spectral clustering. \code{SSC} is 
#' an implementation based on basis pursuit for sparse reconstruction for the 
#' model without systematic noise, which solves
#' \deqn{\textrm{min}_C \|C\|_1\quad\textrm{such that}\quad diag(C)=0,~D=DC}
#' for column-stacked data matrix \eqn{D}. If you are interested in 
#' full implementation of the algorithm with sparse outliers and noise, please 
#' contact the maintainer.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @references 
#' \insertRef{elhamifar_sparse_2009}{T4cluster}
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
#' ## run SSC algorithm with k=2, 3, and 4
#' output2 = SSC(data, k=2)
#' output3 = SSC(data, k=3)
#' output4 = SSC(data, k=4)
#' 
#' ## extract label information
#' lab2 = output2$cluster
#' lab3 = output3$cluster
#' lab4 = output4$cluster
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(3,4))
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=lab2,main="K=2:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=lab2,main="K=2:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=lab2,main="K=2:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=lab2,main="K=2:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=lab3,main="K=3:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=lab3,main="K=3:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=lab3,main="K=3:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=lab3,main="K=3:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=lab4,main="K=4:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=lab4,main="K=4:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=lab4,main="K=4:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=lab4,main="K=4:Axis(2,3)")
#' par(opar)
#' }
#' 
#' @concept subspace
#' @export
SSC <- function(data, k=2){
  ## PREPARE : EXPLICIT INPUTS
  X = prec_input_matrix(data)
  N = base::nrow(X)
  K = max(1, round(k))

  # Step 1. Solve Multiple ADMM Program
  C = array(0,c(N,N))
  for (n in 1:N){
    C[,n] = as.vector(ssc_single_bp(X, n))
  }
  
  # Step 2. Normalize the Columns of C
  for (n in 1:N){
    tgtvec = C[,n]
    C[,n]  = tgtvec/base::max(base::abs(tgtvec))
  }
  
  # Step 3. Construct Similarity Graph
  W = (base::abs(C) + t(base::abs(C)))/2
  
  # Step 4. Spectral Clustering 
  #         Paper says NJW .. okay
  cpprun = sc_normalNJW(W, K, TRUE, 100)

  ## WRAP
  output  = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  # output$eigval  = as.vector(cpprun$values)
  # output$embeds  = cpprun$embeds
  output$algorithm = "SSC"
  return(structure(output, class="T4cluster"))
}

#' @keywords internal
#' @noRd
ssc_single_bp <- function(Y, id){
  n = base::nrow(Y)
  
  myA <- t(Y[-id,])        # p x (n-1)
  myb <- as.vector(Y[id,]) # p
  
  sol <- as.vector(ADMM::admm.bp(myA, myb, maxiter=100)$x)
  if (id==1){
    return(c(0, sol))
  } else if (id>=n){
    return(c(sol, 0))
  } else {
    return(c(sol[1:(id-1)], 0, sol[id:(n-1)]))
  }
}