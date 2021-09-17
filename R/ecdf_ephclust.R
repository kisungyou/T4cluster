#' Hierarchical Agglomerative Clustering for Empirical Distributions
#' 
#' Given \eqn{N} empirical CDFs, perform hierarchical clustering.
#' 
#' @param elist a length-\eqn{N} list of \code{\link[stats]{ecdf}} objects or arrays that can be converted into a numeric vector.
#' @param method agglomeration method to be used. This must be one of \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param ... extra parameters including \describe{
#' \item{type}{(case-insensitive) type of the distance measures (default: \code{"ks"}).}
#' \item{p}{order for the distance for metrics including \code{Wasserstein} and \code{lp} (default: 2).}
#' }
#' 
#' @return an object of \code{hclust} object. See \code{\link[stats]{hclust}} for details. 
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #              3 Types of Univariate Distributions
#' #
#' #    Type 1 : Mixture of 2 Gaussians
#' #    Type 2 : Gamma Distribution
#' #    Type 3 : Mixture of Gaussian and Gamma
#' # -------------------------------------------------------------
#' # generate data
#' myn   = 50
#' elist = list()
#' for (i in 1:10){
#'    elist[[i]] = stats::ecdf(c(rnorm(myn, mean=-2), rnorm(myn, mean=2)))
#' }
#' for (i in 11:20){
#'    elist[[i]] = stats::ecdf(rgamma(2*myn,1))
#' }
#' for (i in 21:30){
#'    elist[[i]] = stats::ecdf(rgamma(myn,1) + rnorm(myn, mean=3))
#' }
#' 
#' # run 'ephclust' with different distance measures
#' eh_ks <- ephclust(elist, type="ks")
#' eh_lp <- ephclust(elist, type="lp")
#' eh_wd <- ephclust(elist, type="wass")
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(eh_ks, main="Kolmogorov-Smirnov")
#' plot(eh_lp, main="L_p")
#' plot(eh_wd, main="Wasserstein")
#' par(opar)
#' }
#' 
#' @concept ecdf
#' @export
ephclust <- function(elist, method=c("single","complete","average","mcquitty",
                                     "ward.D","ward.D2","centroid","median"), ...){
  ###############################################
  # preprocessing
  clist    = elist_converter(elist, fname="ephclust")
  mymethod = match.arg(method)
  
  pars     = list(...)
  pnames   = names(pars)
  if ("type"%in%pnames){
    mytype = match.arg(tolower(pars$type), c("ks","lp","wass"))
  } else {
    mytype = "ks"
  }
  if ("p"%in%pnames){
    myp = as.double(pars$p)
    if (myp < 0){
      stop("* ephclust : 'p' should be a nonnegative number.")
    }
  } else {
    myp = 2.0
  }
  
  ###############################################
  # compute pairwise distance
  dmat = elist_pdist(clist, mytype, myp)
  
  ###############################################
  # apply 'maotai'
  hfun <- utils::getFromNamespace("hidden_hclust","maotai")
  output <- hfun(dmat, mymethod, NULL)
  return(output)
}