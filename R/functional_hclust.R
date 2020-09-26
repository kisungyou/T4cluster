#' Functional Hierarchical Clustering
#' 
#' Given \eqn{N} curves \eqn{\gamma_1 (t), \gamma_2 (t), \ldots, \gamma_N (t) : I \rightarrow \mathbf{R}}, 
#' perform hierarchical agglomerative clustering with \pkg{fastcluster} package's implementation of 
#' the algorithm. Dissimilarity for curves is measured by \eqn{L_p} metric. 
#' 
#' @param fdobj a \code{'fd'} functional data object of \eqn{N} curves by the \pkg{fda} package.
#' @param p an exponent in \eqn{L_p} formalism (default: 2). 
#' @param method agglomeration method to be used. This must be one of \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details. 
#' 
#' @examples
#' # -------------------------------------------------------------
#' #                     two types of curves
#' #
#' # type 1 : sin(x) + perturbation; 20 OF THESE ON [0, 2*PI]
#' # type 2 : cos(x) + perturbation; 20 OF THESE ON [0, 2*PI]
#' # -------------------------------------------------------------
#' ## PREPARE : USE 'fda' PACKAGE
#' #  Generate Raw Data
#' datx = seq(from=0, to=2*pi, length.out=100)
#' daty = array(0,c(100, 40))
#' for (i in 1:20){
#'   daty[,i]    = sin(datx) + rnorm(100, sd=0.1)
#'   daty[,i+20] = cos(datx) + rnorm(100, sd=0.1)
#' }
#' #  Wrap as 'fd' object
#' mybasis <- fda::create.bspline.basis(c(0,2*pi), nbasis=10)
#' myfdobj <- fda::smooth.basis(datx, daty, mybasis)$fd
#' 
#' ## RUN THE ALGORITHM 
#' hcsingle = funhclust(myfdobj, method="single")
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' matplot(datx, daty[,1:20],  type="l", main="Curves Type 1")
#' matplot(datx, daty[,21:40], type="l", main="Curves Type 2")
#' plot(hcsingle, main="hclust with 'single' linkage")
#' par(opar)
#' 
#' @concept functional
#' @export
funhclust <- function(fdobj, p=2, 
                      method = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
                                 "centroid", "median"), members=NULL){
  ## PREPARE
  if (!inherits(fdobj,"fd")){
    stop("* funhclust : input 'fdobj' should be a 'fd' object.")
  }
  myp      = max(1, as.double(p))
  mymethod = match.arg(method)
  
  ## COMPUTE PAIRWISE DISTANCE
  pdmat = stats::as.dist(faux_pdist_lp(fdobj, p=myp))
  
  ## IMPORT THE FUNCTION
  fimport = utils::getFromNamespace("hidden_hclust", "maotai")
  hcout   = fimport(pdmat, mymethod, members)
  return(hcout)
}

# basis <- create.bspline.basis(c(0, 365), nbasis=21, norder=4)
# fdobj <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],basis,
#                       fdnames=list("Day", "Station", "Deg C"))$fd
# fdobj <- smooth.basis(day.5, as.vector(CanadianWeather$dailyAv[,1,"Temperature.C"]),basis)$fd
# 
# data1 <- CanadianWeather$dailyAv[,1,"Temperature.C"]
# data2 <- eval.fd(0:365, fdobj)[,1]
# 
# par(mfrow=c(1,2))
# plot(data1)
# plot(data2)