#' Geodesic Spherical K-Means
#' 
#' Geodesic spherical \eqn{k}-means algorithm is an counterpart of the spherical \eqn{k}-means 
#' algorithm by replacing the cosine similarity with the squared geodesic distance, 
#' which is the great-circle distance under the intrinsic geometry regime 
#' on the unit hypersphere.  If the data is not 
#' normalized, it performs the normalization and proceeds thereafter.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations. If not row-stochastic, each row is normalized to be unit norm.
#' @param k the number of clusters (default: 2).
#' @param ... extra parameters including \describe{
#' \item{init}{initialization method; either \code{"kmeans"} or \code{"gmm"} (default: \code{"kmeans"}).}
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{abstol}{stopping criterion to stop the algorithm (default: \eqn{10^{-8}}).}
#' \item{verbose}{a logical; \code{TRUE} to show iteration history or \code{FALSE} to quiet.}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{cost}{a value of the cost function.}
#' \item{means}{an \eqn{(k\times p)} matrix where each row is a unit-norm class mean. }
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #            clustering with 'household' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' data(household, package="T4cluster")
#' X   = household$data
#' lab = as.integer(household$gender)
#' 
#' ## EXECUTE GSKMEANS WITH VARYING K's
#' vec.rand = rep(0, 9)
#' for (i in 1:9){
#'   clust_i = gskmeans(X, k=(i+1))$cluster
#'   vec.rand[i] = compare.rand(clust_i, lab)
#' }
#' 
#' ## VISUALIZE THE RAND INDEX
#' opar <- par(no.readonly=TRUE)
#' plot(2:10, vec.rand, type="b", pch=19, ylim=c(0.5, 1),
#'      ylab="Rand index",xlab="number of clusters",
#'      main="clustering quality index over varying k's.")
#' par(opar)
#' }
#' 
#' @concept sphere
#' @export
gskmeans <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUTS
  mydata = prec_input_sphere(as.matrix(data))
  myk    = max(1, round(k))
  
  ## PREPARE : IMPLICIT ONES
  params  = list(...)
  pnames  = names(params)
  
  if ("maxiter"%in%pnames){    myiter = max(5, round(params$maxiter))  } else {   myiter = 10  }
  if ("abstol"%in%pnames){    myeps = max(params$abstol, .Machine$double.eps)  } else {    myeps = sqrt(.Machine$double.eps)  }
  if ("init"%in% pnames){
    myinit = match.arg(tolower(params$init), c("kmeans","gmm"))
  } else {
    myinit = "kmeans"
  }
  if ("verbose"%in%pnames){
    myprint = as.logical(params$verbose)
  } else {
    myprint = FALSE
  }
  
  ## RUN
  cpprun = sp_gskmeans(mydata, myk, myinit, myiter, myeps, myprint)
  
  ## WRAP
  output  = list()
  output$cluster   = round(as.vector(cpprun$cluster+1))
  output$cost      = as.double(cpprun$cost)
  output$means     = cpprun$means
  output$algorithm = "gskmeans"
  return(structure(output, class="T4cluster"))   
}

#' @return a list containing\describe{
#' \item{data}{an \eqn{(n\times 2)} data matrix.}
#' \item{label}{a length-\eqn{n} vector(factor) for class labels.}
#' } 

