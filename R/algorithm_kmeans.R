#' K-Means Clustering
#' 
#' \eqn{K}-means algorithm we provide is a wrapper to the \pkg{Armadillo}'s k-means routine.
#' Two types of initialization schemes are employed. Please see the parameters section for more details.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param ... extra parameters including \describe{
#' \item{init}{initialization method; either \code{"random"} for random initialization, or \code{"plus"} for k-means++ starting.}
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{nstart}{the number of random initializations (default: 5).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{mean}{a \eqn{(k\times p)} matrix where each row is a class mean.}
#' \item{wcss}{within-cluster sum of squares (WCSS).}
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
#' cl2 = kmeans(X, k=2)$cluster
#' cl3 = kmeans(X, k=3)$cluster
#' cl4 = kmeans(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="k-means: k=2")
#' plot(X2d, col=cl3, pch=19, main="k-means: k=3")
#' plot(X2d, col=cl4, pch=19, main="k-means: k=4")
#' par(opar)
#' 
#' @references 
#' \insertRef{sanderson_armadillo_2016c}{T4cluster}
#' 
#' @concept algorithm
#' @export
kmeans <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUTS
  mydata = prec_input_matrix(data)
  myk    = max(1, round(k))
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  myiter   = ifelse(("maxiter" %in% pnames), max(10, round(pars$maxiter)), 10)
  myinit   = ifelse(("init"%in%pnames), match.arg(tolower(pars$init),c("random","plus")),"plus")
  mynstart = ifelse(("nstart"%in%pnames), max(5,round(pars$nstart)), 5)
  
  ## MULTIPLE STARTS
  rec_class = list()
  rec_kmeans = list()
  if (all(myinit=="random")){
    for (i in 1:mynstart){
      rec_kmeans[[i]] = arma_kmeans_random(t(data), myk, myiter)
    }
  } else {
    dmatrix = cpp_pdistMP(mydata, 2, 4)
    distobj = stats::as.dist(dmatrix)
    for (i in 1:mynstart){
      tmplab = extra_kmeanspp(distobj, k=myk)$id.center
      rec_kmeans[[i]] = arma_kmeans_kmeanspp(t(data), t(data[tmplab,]), myk, myiter)
    }
  }
  
  ## FIND THE BEST ONE
  vec_SSE = rep(0,mynstart)
  for (i in 1:mynstart){
    tgtrun     = rec_kmeans[[i]]
    tmplabel   = apply(tgtrun$pdmat, 1, which.min)
    vec_SSE[i] = kmeans_SSE(tgtrun$pdmat, tmplabel)
  }
  optrun = rec_kmeans[[which.min(vec_SSE)]]
  
  ## SELECT, WRAP, AND RUN
  output = list()
  output$cluster   = base::apply(optrun$pdmat, 1, which.min)
  output$means     = optrun$means
  output$wcss      = min(vec_SSE)
  output$algorithm ="kmeans"
  return(structure(output, class="T4cluster"))
}


# extra for kmeans --------------------------------------------------------
#' @keywords internal
#' @noRd
kmeans_SSE <- function(pdmat, label){
  n = nrow(pdmat)
  k = ncol(pdmat)
  output = 0
  for (i in 1:n){
    output = output + (as.double(pdmat[i,label[i]])^2)
  }
  return(output)
}
