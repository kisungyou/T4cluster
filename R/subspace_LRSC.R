#' Low-Rank Sparse Clustering
#' 
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
#' ## run LRSC algorithm with k=2,3,4 with relaxed/exact solvers
#' out2rel = LRSC(data, k=2, type="relaxed")
#' out3rel = LRSC(data, k=3, type="relaxed")
#' out4rel = LRSC(data, k=4, type="relaxed")
#' 
#' out2exc = LRSC(data, k=2, type="exact")
#' out3exc = LRSC(data, k=3, type="exact")
#' out4exc = LRSC(data, k=4, type="exact")
#' 
#' ## extract label information
#' lab2rel = out2rel$cluster
#' lab3rel = out3rel$cluster
#' lab4rel = out4rel$cluster
#' 
#' lab2exc = out2exc$cluster
#' lab3exc = out3exc$cluster
#' lab4exc = out4exc$cluster
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3))
#' plot(dat2, pch=19, cex=0.9, col=lab2rel, main="LRSC Relaxed:K=2")
#' plot(dat2, pch=19, cex=0.9, col=lab3rel, main="LRSC Relaxed:K=3")
#' plot(dat2, pch=19, cex=0.9, col=lab4rel, main="LRSC Relaxed:K=4")
#' plot(dat2, pch=19, cex=0.9, col=lab2exc, main="LRSC Exact:K=2")
#' plot(dat2, pch=19, cex=0.9, col=lab3exc, main="LRSC Exact:K=3")
#' plot(dat2, pch=19, cex=0.9, col=lab4exc, main="LRSC Exact:K=4")
#' par(opar)
#' }
#' 
#' 
#' @concept subspace
#' @export
LRSC <- function(data, k=2, type=c("relaxed","exact"), tau=1.0){
  ## PREPARE : EXPLICIT INPUTS
  X = prec_input_matrix(data)
  myk     = max(1, round(k))
  mytau   = max(sqrt(.Machine$double.eps), as.double(tau))
  myexact = match.arg(type)
  
  ## RUN EVERYTHING IN C++
  cpprun = cpp_LRSC(X, myk, myexact, mytau)

  ## WRAP
  output  = list()
  output$cluster = round(as.vector(cpprun$labels+1))
  output$algorithm = "LRSC"
  return(structure(output, class="T4cluster"))
}