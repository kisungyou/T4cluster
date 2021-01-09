#' Weighted GMM by Gebru et al. (2016)
#' 
#' When each observation \eqn{x_i} is associated with a weight \eqn{w_i > 0}, 
#' modifying the GMM formulation is required. Gebru et al. (2016) proposed a method 
#' to use scaled covariance based on an observation that
#' \deqn{\mathcal{N}\left(x\vert \mu, \Sigma\right)^w \propto \mathcal{N}\left(x\vert \mu, \frac{\Sigma}{w}\right)} 
#' by considering the positive weight as a role of precision. Currently, 
#' we provide a method with fixed weight case only while the paper also considers 
#' a Bayesian formalism on the weight using Gamma distribution.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param weight a positive weight vector of length \eqn{n}. If \code{NULL} (default), uniform weight is set.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 10).}
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
#' ## CLUSTERING WITH DIFFERENT K VALUES
#' cl2 = gmm16G(X, k=2)$cluster
#' cl3 = gmm16G(X, k=3)$cluster
#' cl4 = gmm16G(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="gmm16G: k=2")
#' plot(X2d, col=cl3, pch=19, main="gmm16G: k=3")
#' plot(X2d, col=cl4, pch=19, main="gmm16G: k=4")
#' par(opar)
#' 
#' @references 
#' \insertRef{gebru_em_2016}{T4cluster}
#' 
#' @concept algorithm
#' @export
gmm16G <- function(data, k=2, weight=NULL, ...){ # currently, fixed method is only provided.
  ## PREPARE : EXPLICIT INPUT
  mydata  = prec_input_matrix(data)
  myk     = max(1, round(k))
  myn     = base::nrow(mydata)
  if ((length(weight)==0)&&(is.null(weight))){
    myweight = rep(1/myn, myn)
  } else {
    myweight = as.vector(weight)
  }
  if ((any(myweight<=0))||(length(myweight)!=myn)){
    stop(paste0("* gmm16G : weight vector should be of length ",myn," of positive values."))
  }
  myweight = (myweight/base::sum(myweight))*myn # this makes normalization ----------
  
  ## PREPARE : IMPLICIT
  params  = list(...)
  pnames  = names(params)
  
  myiter   = round(max(5, ifelse(("maxiter"%in%pnames), params$maxiter, 10)))
  mydiag   = as.logical(ifelse(("usediag"%in%pnames), params$usediag, FALSE))

  ## RUN
  cpprun = gmm_16Gfix(data, myk, myweight, myiter, mydiag)
  
  ## WRAP AND RETURN
  output = list()
  output$cluster  = as.vector(cpprun$cluster+1)
  output$mean     = cpprun$means
  output$variance = cpprun$covs
  output$weight   = as.vector(cpprun$weight)
  output$loglkd   = as.double(cpprun$loglkd)
  output$algorithm = "gmm16G"
  return(structure(output, class="T4cluster"))  
}
