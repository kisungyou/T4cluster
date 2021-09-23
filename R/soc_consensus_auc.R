#' Consensus Determination of K : Original Method
#' 
#' @param clist a list of consensus matrices wrapped as \code{T4cluster:consensus} objects.
#' 
#' @return a \code{data.frame} consisting of \describe{
#' \item{\code{nclust}}{numbers of clusters.}
#' \item{\code{delta}}{delta values as noted in the paper.}
#' }
#' 
#' @examples 
#' \donttest{
#' 
#' }
#' 
#' 
#' 
#' @concept soc
#' @export
consensus.auc <- function(clist){
  ## PREPARE
  check_clist(clist, "consensus.auc")
  my_K   = 2:(length(clist)+1)

  ## COMPUTE
  nobj   = length(clist)
  my_auc = rep(0, nobj)
  for (i in 1:nobj){
    tgt_mat   = clist[[i]]
    tgt_ecdf  = stats::ecdf(as.vector(tgt_mat[upper.tri(tgt_mat)]))
    my_auc[i] = ecdf_auc(tgt_ecdf)
  }
  
  ## COMPUTE DELTA
  vec_delta = rep(0, nobj)
  vec_delta[1] = my_auc[1]
  for (i in 2:nobj){
    vec_delta[i] = (my_auc[i] - my_auc[i-1])/my_auc[i-1]
  }
  
  ## WRAP FOR RETURN
  output = list()
  output$record = data.frame(nclust=my_K, delta=vec_delta)
  return(output)
}


# X = T4cluster::genSMILEY(n=100)
# data = X$data
# labs = X$label
# 
# fun_phate <- utils::getFromNamespace("hidden_PHATE","maotai")
# run_phate <- fun_phate(dist(data), nbdk=5, alpha=2)
# pseudoX   <- base::sqrt(run_phate$P)
# 
# runs = list()
# nmax = 10
# for (i in 1:nmax){
#   runs[[i]] = ccphate(pseudoX, k=(i+1))
#   print(paste0("iteration ",i," complete...."))
# }
# 
# cclist = list()
# for (i in 1:nmax){
#   cclist[[i]] = runs[[i]]$consensus
# }
# cobj = consensus.auc(cclist)$record
# plot(cobj$nclust, cobj$delta, "b", xaxt='n')
# axis(1, at=cobj$nclust)
