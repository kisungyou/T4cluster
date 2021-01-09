#' Regularized GMM by Ruan et al. (2011)
#' 
#' Ruan et al. (2011) proposed a regularized covariance estimation by 
#' graphical lasso to cope with high-dimensional scenario where conventional 
#' GMM might incur singular covariance components. Authors proposed to use 
#' \eqn{\lambda} as a regularization parameter as normally used in 
#' sparse covariance/precision estimation problems and suggested to use the 
#' model with the smallest BIC values. 
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param lambda regularization parameter for graphical lasso (default: 1).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{nstart}{the number of random initializations (default: 5).}
#' \item{usediag}{a logical; covariances are diagonal if \code{TRUE}, or full covariances are returned for \code{FALSE} (default: \code{FALSE}).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).}
#' \item{mean}{a \eqn{(k\times p)} matrix where each row is a class mean.}
#' \item{variance}{a \eqn{(p\times p\times k)} array where each slice is a class covariance.}
#' \item{weight}{a length-\eqn{k} vector of class weights that sum to 1.}
#' \item{loglkd}{log-likelihood of the data for the fitted model.}
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
#' ## COMPARE WITH STANDARD GMM
#' cl.gmm = gmm(X, k=3)$cluster
#' cl.11Rf = gmm11R(X, k=3)$cluster
#' cl.11Rd = gmm11R(X, k=3, usediag=TRUE)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(X2d, col=cl.gmm,  pch=19, main="standard GMM")
#' plot(X2d, col=cl.11Rf, pch=19, main="gmm11R: full covs")
#' plot(X2d, col=cl.11Rd, pch=19, main="gmm11R: diagonal covs")
#' par(opar)
#' 
#' @references 
#' \insertRef{ruan_regularized_2011}{T4cluster}
#' 
#' @concept algorithm
#' @export
gmm11R <- function(data, k=2, lambda=1, ...){
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_matrix(data)
  myk     = max(1, round(k))
  mylbd   = max(0.1, as.double(lambda))
  
  ## PREPARE : IMPLICIT
  params  = list(...)
  pnames  = names(params)
  
  myiter   = round(max(5, ifelse(("maxiter"%in%pnames), params$maxiter, 10)))
  mydiag   = as.logical(ifelse(("usediag"%in%pnames), params$usediag, FALSE))
  mynstart = max(1,round(ifelse(("nstart"%in%pnames), params$nstart, 5)))
  
  ## MULTIPLE RUNS
  runs.list = list()
  runs.lkd  = rep(0,mynstart)
  for (i in 1:mynstart){
    runs.list[[i]] = gmm_11R(mydata, myk, mylbd, myiter, mydiag)
    runs.lkd[i]    = runs.list[[i]]$loglkd
  }
  cpprun = runs.list[[which.max(runs.lkd)]]
  
  ## WRAP AND RETURN
  output = list()
  output$cluster  = as.vector(cpprun$cluster+1)
  output$mean     = cpprun$means
  output$variance = cpprun$covs
  output$weight   = as.vector(cpprun$weight)
  output$loglkd   = as.double(cpprun$loglkd)
  output$algorithm = "gmm11R"
  return(structure(output, class="T4cluster"))  
}


  