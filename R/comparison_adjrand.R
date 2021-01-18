#' (+) Adjusted Rand Index
#' 
#' Compute Adjusted Rand index between two clusterings. Please note that the 
#' value can yield negative value.
#' 
#' @seealso \code{\link{compare.rand}}
#' 
#' @param x 1st cluster label vector of length-\eqn{n}.
#' @param y 2nd cluster label vector of length-\eqn{n}.
#' 
#' @return Adjusted Rand Index value.
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #         true label vs. clustering with 'iris' dataset 
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(iris)
#' X   = as.matrix(iris[,1:4])
#' lab = as.integer(as.factor(iris[,5]))
#' 
#' ## CLUSTERING WITH DIFFERENT K VALUES
#' vec_k  = 2:7
#' vec_cl = list()
#' for (i in 1:length(vec_k)){
#'   vec_cl[[i]] = T4cluster::kmeans(X, k=round(vec_k[i]))$cluster
#' }
#' 
#' ## COMPUTE COMPARISON INDICES
#' vec_comp = rep(0, length(vec_k))
#' for (i in 1:length(vec_k)){
#'   vec_comp[i] = compare.adjrand(vec_cl[[i]], lab)
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(vec_k, vec_comp, type="b", lty=2, xlab="number of clusters", 
#'      ylab="comparison index", main="Adjusted Rand Index with true k=3")
#' abline(v=3, lwd=2, col="red")
#' par(opar)
#' }
#' 
#' @concept comparison
#' @export
compare.adjrand <- function(x, y){
  # Preprocessing
  fname = "compare.adjrand"
  aname = "adjrand"
  
  x = as.integer(as.factor(x))
  y = as.integer(as.factor(y))
  prec_twolabel(x, y, fname)
  
  # Compute
  return(as.double(mclustcomp::mclustcomp(x, y, aname)$scores[1]))
}