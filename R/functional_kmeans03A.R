#' Functional K-Means Clustering by Abraham et al. (2003)
#' 
#' Given \eqn{N} curves \eqn{\gamma_1 (t), \gamma_2 (t), \ldots, \gamma_N (t) : I \rightarrow \mathbf{R}}, 
#' perform \eqn{k}-means clustering on the coefficients from the functional data expanded by 
#' B-spline basis. Note that in the original paper, authors used B-splines as the choice of basis 
#' due to nice properties. However, we allow other types of basis as well for convenience.
#' 
#' @param fdobj a \code{'fd'} functional data object of \eqn{N} curves by the \pkg{fda} package.
#' @param k the number of clusters (default: 2).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{nstart}{the number of random initializations (default: 5).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{N} vector of class labels (from \eqn{1:k}).} 
#' \item{mean}{a \code{'fd'} object of \eqn{k} mean curves.}
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @examples
#' # -------------------------------------------------------------
#' #                     two types of curves
#' #
#' # type 1 : sin(x) + perturbation; 20 OF THESE ON [0, 2*PI]
#' # type 2 : cos(x) + perturbation; 20 OF THESE ON [0, 2*PI]
#' # type 3 : sin(x) + cos(0.5x)   ; 20 OF THESE ON [0, 2*PI]
#' # -------------------------------------------------------------
#' ## PREPARE : USE 'fda' PACKAGE
#' #  Generate Raw Data
#' datx = seq(from=0, to=2*pi, length.out=100)
#' daty = array(0,c(100, 60))
#' for (i in 1:20){
#'   daty[,i]    = sin(datx) + rnorm(100, sd=0.5)
#'   daty[,i+20] = cos(datx) + rnorm(100, sd=0.5)
#'   daty[,i+40] = sin(datx) + cos(0.5*datx) + rnorm(100, sd=0.5)
#' }
#' #  Wrap as 'fd' object
#' mybasis <- fda::create.bspline.basis(c(0,2*pi), nbasis=10)
#' myfdobj <- fda::smooth.basis(datx, daty, mybasis)$fd
#' 
#' ## RUN THE ALGORITHM WITH K=2,3,4
#' fk2 = funkmeans03A(myfdobj, k=2)
#' fk3 = funkmeans03A(myfdobj, k=3)
#' fk4 = funkmeans03A(myfdobj, k=4)
#' 
#' ## FUNCTIONAL PCA FOR VISUALIZATION
#' embed = fda::pca.fd(myfdobj, nharm=2)$score
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(embed, col=fk2$cluster, pch=19, main="K=2")
#' plot(embed, col=fk3$cluster, pch=19, main="K=3")
#' plot(embed, col=fk4$cluster, pch=19, main="K=4")
#' par(opar)
#' 
#' @references 
#' \insertRef{abraham_unsupervised_2003}{T4cluster}
#' 
#' @concept functional
#' @export
funkmeans03A <- function(fdobj, k=2, ...){
  ## PREPARE : EXPLICIT INPUTS
  if (!inherits(fdobj,"fd")){
    stop("* funkmeans03A : input 'fdobj' should be a 'fd' object.")
  }
  myk = max(1, round(k))
  if (!(fdobj$basis$type=="bspline")){
    warning("* funkmeans03A : original paper suggested to use B-splines, but we proceed anyway.")
  }
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  myiter   = ifelse(("maxiter" %in% pnames), max(10, round(pars$maxiter)), 10)
  mynstart = ifelse(("nstart"%in%pnames), max(5,round(pars$nstart)), 5)
  
  ## MULTIPLE STARTS
  rec_kmeans = list()
  for (i in 1:mynstart){
    rec_kmeans[[i]] = arma_kmeans_random(as.matrix(fdobj$coefs), myk, myiter)
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
  meanobj = fdobj
  meanobj$coefs = t(optrun$means)
  rownames(meanobj$coefs) = rownames(fdobj$coefs)
  
  output = list()
  output$cluster   = base::apply(optrun$pdmat, 1, which.min)
  output$means     = meanobj
  output$algorithm ="funkmeans03A"
  return(structure(output, class="T4cluster"))
}