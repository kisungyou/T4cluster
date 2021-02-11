#' Ensembles of K-Subspaces
#' 
#' Ensembles of K-Subspaces method exploits multiple runs of K-Subspace Clustering and 
#' uses consensus framework to aggregate multiple clustering results 
#' to mitigate the effect of random initializations. When the results are merged, 
#' it zeros out \eqn{n-q} number of values in a co-occurrence matrix. The paper 
#' suggests to use large number of runs (\code{B}) where each run may not require 
#' large number of iterations (\code{iter}) since the main assumption of the 
#' algorithm is to utilize multiple partially-correct information. At the extreme case, 
#' iteration \code{iter} may be set to 0 for which the paper denotes it as EKSS-0.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param d candidate dimension for each subspace (default: 2).
#' @param q threshold; the number of smaller values to be zeroed out (default: 0.75*\eqn{n}).
#' @param B the number of ensembles/runs (default: 500).
#' @param iter the number of iteration for each run (default: 0).
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
#' ## run EKSS algorithm with k=2,3,4 with EKSS-0 and 5 iterations
#' out2zero = EKSS(data, k=2)
#' out3zero = EKSS(data, k=3)
#' out4zero = EKSS(data, k=4)
#' 
#' out2iter = EKSS(data, k=2, iter=5)
#' out3iter = EKSS(data, k=3, iter=5)
#' out4iter = EKSS(data, k=4, iter=5)
#' 
#' ## extract label information
#' lab2zero = out2zero$cluster
#' lab3zero = out3zero$cluster
#' lab4zero = out4zero$cluster
#' 
#' lab2iter = out2iter$cluster
#' lab3iter = out3iter$cluster
#' lab4iter = out4iter$cluster
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3))
#' plot(dat2, pch=19, cex=0.9, col=lab2zero, main="EKSS-0:K=2")
#' plot(dat2, pch=19, cex=0.9, col=lab3zero, main="EKSS-0:K=3")
#' plot(dat2, pch=19, cex=0.9, col=lab4zero, main="EKSS-0:K=4")
#' plot(dat2, pch=19, cex=0.9, col=lab2iter, main="EKSS iter:K=2")
#' plot(dat2, pch=19, cex=0.9, col=lab3iter, main="EKSS iter:K=3")
#' plot(dat2, pch=19, cex=0.9, col=lab4iter, main="EKSS iter:K=4")
#' par(opar)
#' }
#' 
#' @references 
#' Lipor J, Hong D, Tan YS, Balzano L (2021). “Subspace Clustering Using Ensembles of \eqn{K}-Subspaces.” arXiv:1709.04744.
#' 
#' @concept subspace
#' @export
EKSS <- function(data, k=2, d=2, q=floor(nrow(data)*0.75), B=500, iter=0){
  ## PREPARE : EXPLICIT INPUTS
  X = prec_input_matrix(data)
  K = max(1, round(k))
  N = base::nrow(X)
  mydim  = max(1, round(d))
  myq    = round(q)
  myB    = max(10, round(B))
  myiter = max(0, round(iter))
  
  ## RUN MULTIPLE TIMES
  ## not parallel yet.. but will.. someday.. but we will stack as columns
  labels = array(0,c(N,myB))
  for (it in 1:myB){
    labels[,it] = EKSS_run_main(X, K, mydim, myiter)
  }
  
  ## CONSTRUCT AFFINITY MATRIX & THRESHOLD
  Anaive = cpp_EKSS_affinity(labels)
  Athr   = EKSS_thresh(Anaive, myq)
  
  ## APPLY SPECTRAL CLUSTERING OF NJW
  cpp_njw = sc_normalNJW(Athr, K, TRUE, 100)
  
  ## WRAP
  output  = list()
  output$cluster = round(as.vector(cpp_njw$labels+1))
  output$algorithm = "EKSS"
  return(structure(output, class="T4cluster"))
}



# auxiliary functions for EKSS --------------------------------------------
#' @keywords internal
#' @noRd
EKSS_thresh <- function(A, q){
  N  = base::nrow(A)
  Nq = round(N-q)
  
  Zrow = array(0,c(N,N))
  Zcol = array(0,c(N,N))
  
  for (i in 1:N){ # process rows
    tgt = as.vector(A[i,])
    tgt[order(tgt)[1:Nq]] = 0
    Zrow[i,] = tgt
  }
  for (j in 1:N){
    tgt = as.vector(A[,j])
    tgt[order(tgt)[1:Nq]] = 0
    Zcol[,j] = tgt
  }
  
  Abar = (Zrow+Zcol)/2
  return(Abar)
}
#' @keywords internal
#' @noRd
EKSS_run_main <- function(X, K, d, iter){
  if (iter < 1){
    cpprun = cpp_EKSS_0(X,K,d)  
  } else {
    cpprun = cpp_EKSS_T(X,K,d,iter)
  }
  output = as.integer(as.factor(as.vector(cpprun)))
  return(output)
}