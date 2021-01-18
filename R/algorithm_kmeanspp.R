#' K-Means++ Clustering
#' 
#' \eqn{K}-means++ algorithm is usually used as a fast initialization scheme, though 
#' it can still be used as a standalone clustering algorithms by first choosing the 
#' centroids and assign points to the nearest centroids.
#' 
#' @param data an \eqn{(n \times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).}
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
#' cl2 = kmeanspp(X, k=2)$cluster
#' cl3 = kmeanspp(X, k=3)$cluster
#' cl4 = kmeanspp(X, k=4)$cluster
#' 
#' ## VISUALIZATION
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X2d, col=lab, pch=19, main="true label")
#' plot(X2d, col=cl2, pch=19, main="k-means++: k=2")
#' plot(X2d, col=cl3, pch=19, main="k-means++: k=3")
#' plot(X2d, col=cl4, pch=19, main="k-means++: k=4")
#' par(opar)
#' 
#' @references 
#' \insertRef{arthur_k-means++:_2007}{T4cluster}
#' 
#' @concept algorithm
#' @export
kmeanspp <- function(data, k=2){
  ## PREPARE : EXPLICIT INPUTS
  mydata = prec_input_matrix(data)
  myk    = max(1, round(k))
  
  ## COMPUTE THE LABEL
  tmpout = extra_kmeanspp(mydata, k=myk)
  
  ## WRAP AND RETURN
  output = list()
  output$cluster = tmpout$cluster
  output$algorithm = "kmeanspp"
  return(structure(output, class="T4cluster"))
  return(output)
}


# extra functions ---------------------------------------------------------
#' @keywords internal
#' @noRd
extra_kmeanspp <- function(x, k=2){ # covers both case & return center id's.
  # prepare
  if (is.matrix(x)){
    x = cpp_pdistMP(x, 2, 4)
  } else {
    x = as.matrix(x)
  }
  n  = nrow(x)
  K  = round(k)
  if (K >= n){
    stop("* kmeanspp : 'k' should be smaller than the number of observations.")
  }
  if (K < 2){
    stop("* kmeanspp : 'k' should be larger than 1.")
  }
  id.now = 1:n
  
  # Computation
  #   initialize
  id.center = base::sample(id.now, 1)
  id.now    = base::setdiff(id.now, id.center)
  #   iterate
  for (i in 1:(K-1)){
    # compute distance to the nearest
    tmpdmat = x[id.now, id.center]
    if (i==1){
      d2vec = as.vector(tmpdmat)^2
      d2vec = d2vec/base::sum(d2vec)
    } else {
      d2vec = as.vector(base::apply(tmpdmat, 1, base::min))^2
      d2vec = d2vec/base::sum(d2vec)
    }
    # sample one
    id.tmp = base::sample(id.now, 1, prob=d2vec)
    # update
    id.center = c(id.center, id.tmp)
    id.now    = base::setdiff(id.now, id.tmp)
  }
  #   let's compute label
  dmat    = x[,id.center]
  cluster = base::apply(dmat, 1, base::which.min)
  
  ##################################################
  # Return
  output = list()
  output$cluster   = cluster
  output$id.center = id.center
  return(output)
}

  