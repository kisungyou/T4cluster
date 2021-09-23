#' Consensus Determination of K : PAC method
#' 
#' 
#' 
#' @concept soc
#' @export
consensus.pac <- function(clist, lb=0.1, ub=0.9){
  ## PREPARE
  check_clist(clist, "consensus.pac")
  my_K   = 2:(length(clist)+1)
  my_eps = sqrt(.Machine$double.eps)
  my_lb  = max(my_eps, as.double(lb))
  my_ub  = min(as.double(ub), 1-my_eps)
  
  ## COMPUTE
  nobj   = length(clist)
  my_pac = rep(0, nobj)
  my_n   = base::nrow(clist[[1]])
  for (i in 1:nobj){
    tgt_mat   = clist[[i]]
    tgt_tri   = as.vector(tgt_mat[upper.tri(tgt_mat)])
    tgt_ecdf  = stats::ecdf(tgt_tri)
    my_pac[i] = tgt_ecdf(my_ub)-tgt_ecdf(my_lb)
    
    
    # rPAC 
    # for(i in 1: length(xvec)) {
    #   M <- ml[[i]]
    #   Fn <- ecdf(M[lower.tri(M)])
    #   # calculate PAC
    #   PAC[i] <- Fn(x2) - Fn(x1)
    #   # calculate proportion of zeroes in M
    #   prop_zeroes[i] <- sum(M==0)/length(M)
    #   # calculate relative PAC
    #   rPAC[i] <- PAC[i]/(1-prop_zeroes[i])
    # }
    # prop_zero = sum(tgt_tri<.Machine$double.eps)/length(tgt_tri)
    # my_pac[i] =  (tgt_ecdf(my_rb)-tgt_ecdf(my_lb))/(1-prop_zero)
  }
  
  ## WRAP FOR RETURN
  output        = list()
  output$K      = my_K[which.min(my_pac)]
  output$record = data.frame(nclust=my_K, pac=my_pac)
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
# x11()
# cobj = consensus.pac(cclist, lb=0.1, ub=0.9)$record
# plot(cobj$nclust, cobj$pac, "b", xaxt='n')
# axis(1, at=cobj$nclust)