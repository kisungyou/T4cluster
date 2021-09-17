#' K-Medoids Clustering for Empirical Distributions
#' 
#' 
#' 
#' @concept ecdf
#' @export
epmedoids <- function(elist, k=2, ...){
  ## PREPROCESSING
  clist    = elist_converter(elist, fname="epmedoids")
  myk      = max(1, round(k))
  
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
      stop("* epmedoids : 'p' should be a nonnegative number.")
    }
  } else {
    myp = round(2)
  }
  
  ## COMPUTE PAIRWISE DISTANCES
  dmat = elist_pdist(clist, mytype, myp)
  
  ## RUN K-MEDOIDS
  func.import  = utils::getFromNamespace("hidden_kmedoids", "maotai")
  obj.kmedoids = func.import(dmat, nclust=myk) 
  
  ## WRAP AND RETURN
  output = list()
  output$medoids = obj.kmedoids$id.med
  output$cluster = as.integer(obj.kmedoids$clustering)
  output$algorithm ="epmedoids"
  return(structure(output, class="T4cluster"))
}