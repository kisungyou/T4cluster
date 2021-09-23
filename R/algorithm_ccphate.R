#' Consensus Clustering with PHATE Geometry
#' 
#' 
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of clusters (default: 2).
#' @param ... extra parameters including \describe{
#' \item{nnbd}{size of nearest neighborhood (default: 5).}
#' \item{alpha}{decay parameter for Gaussian kernel exponent (default: 2).}
#' \item{ratio}{percentage of samples to be used for each resampling run with the minimum is set as 0.25 (default: 0.9).}
#' \item{niter}{number of resampling runs where the minimum is set 50 (default: 100).}
#' }
#' 
#' @concept algorithm
#' @export
ccphate <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUTS
  mydata = prec_input_matrix(data)
  myk    = max(1, round(k))
  myn    = base::nrow(data)
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  if ("nnbd"%in%pnames){
    my_nnbd = max(1, round(pars$nnbd))
  } else {
    my_nnbd = 5
  }
  if ("alpha"%in%pnames){
    my_alpha = max(1, as.double(pars$alpha))
  } else {
    my_alpha = 2
  }
  if ("ratio"%in%pnames){
    my_ratio = min(max(0.25, as.double(pars$ratio)), 0.99999999)
  } else {
    my_ratio = 0.90
  }
  my_subsize = round(myn*my_ratio)
  if ("niter" %in% pnames){
    my_niter = max(50, round(pars$niter))
  } else {
    my_niter = 100
  }
  
  ## RUN PHATE FOR OPTIMAL ROW-STOCHASTIC MATRIX RETRIEVAL
  if ((is.matrix(mydata))&&(all(rowSums(mydata)==1))&&(nrow(mydata)==ncol(mydata))){
    pseudoX <- base::sqrt(mydata)
  } else {
    distX     <- stats::dist(mydata)
    fun_phate <- utils::getFromNamespace("hidden_PHATE","maotai")
    run_phate <- fun_phate(distX, nbdk=my_nnbd, alpha=my_alpha)
    pseudoX   <- base::sqrt(run_phate$P)
  }
  
  ## RUN MULTIPLE SPHERICAL K-MEANS WITH EXTRINSIC DISTANCE
  myinit  = "kmeans"
  myprint = FALSE

  M_h = array(0,c(myn,myn))
  I_h = array(0,c(myn,myn))
  
  for (it in 1:my_niter){
    # subsampling
    select_id  = base::sample(1:myn, my_subsize, replace=FALSE)
    select_dat = pseudoX[select_id,]
    select_cpp = sp_spkmeans(select_dat, myk, myinit, 100, 1e-8, myprint) 
    select_lab = round(as.vector(select_cpp$cluster)+1)
    
    # update I_h
    for (i in 1:(my_subsize-1)){
      id_i = select_id[i]
      for (j in (i+1):my_subsize){
        id_j = select_id[j]
        I_h[id_i,id_j] = I_h[id_i,id_j] + 1
        I_h[id_j,id_i] = I_h[id_i,id_j]
      }
    }
    # update M_h
    for (i in 1:(my_subsize-1)){
      id_i = select_id[i]
      for (j in (i+1):my_subsize){
        id_j = select_id[j]
        if (select_lab[i]==select_lab[j]){
          M_h[id_i, id_j] = M_h[id_i, id_j] + 1
          M_h[id_j, id_i] = M_h[id_i, id_j]
        }
      }
    }
  }
  
  ## CONSENSUS MATRIX
  M = M_h/I_h
  diag(M) = 1.0
  M[is.nan(M)] = 0.0
  
  ## APPLY K-MEDOIDS FOR A SINGLE CLUSTER OUTPUT : Change to Spectral Clustering
  pseudoD = stats::as.dist(1-M)
  medoid_fun = utils::getFromNamespace("hidden_kmedoids", "maotai")
  medoid_run = medoid_fun(pseudoD, myk)$clustering
  recovered  = as.vector(medoid_run)
  # recovered = as.vector(sc_normalSM(M, myk, TRUE, 100)$labels)+1
  
  ## WRAP
  output  = list()
  output$cluster   = recovered
  output$consensus = structure(M, class="T4cluster:consensus")
  output$algorithm = "ccphate"
  return(structure(output, class="T4cluster"))
}

# graphics.off()
# X = T4cluster::genSMILEY(n=100)
# data = X$data
# labs = X$label
# 
# fun_phate <- utils::getFromNamespace("hidden_PHATE","maotai")
# run_phate <- fun_phate(dist(data), nbdk=5, alpha=2)
# pseudoX   <- base::sqrt(run_phate$P)
# 
# run3 = ccphate(pseudoX, k=3)
# run4 = ccphate(pseudoX, k=4)
# run5 = ccphate(pseudoX, k=5)
# run6 = ccphate(pseudoX, k=6)
# 
# x11()
# par(mfrow=c(4,2), pty="s")
# image(run3$consensus, main="K=3"); plot(data, col=run3$cluster, pch=19)
# image(run4$consensus, main="K=4"); plot(data, col=run4$cluster, pch=19)
# image(run5$consensus, main="K=5"); plot(data, col=run5$cluster, pch=19)
# image(run6$consensus, main="K=6"); plot(data, col=run6$cluster, pch=19)
