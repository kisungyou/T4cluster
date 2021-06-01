#' Spherical K-Means Clustering 
#' 
#' Spherical \eqn{k}-means algorithm performs clustering for the data residing 
#' on the unit hypersphere with the cosine similarity. If the data is not 
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
#' ## EXECUTE SPKMEANS WITH VARYING K's
#' vec.rand = rep(0, 9)
#' for (i in 1:9){
#'   clust_i = spkmeans(X, k=(i+1))$cluster
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
#' @references 
#' I. S. Dhillon and D. S. Modha (2001). "Concept decompositions for large sparse text data using clustering." \emph{Machine Learning}, \strong{42}:143–175.
#' 
#' 
#' @concept sphere
#' @export
spkmeans <- function(data, k=2, ...){
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
  cpprun = sp_spkmeans(mydata, myk, myinit, myiter, myeps, myprint)
  
  ## WRAP
  output  = list()
  output$cluster   = round(as.vector(cpprun$cluster+1))
  output$cost      = as.double(cpprun$cost)
  output$means     = cpprun$means
  output$algorithm = "spkmeans"
  return(structure(output, class="T4cluster"))   
}