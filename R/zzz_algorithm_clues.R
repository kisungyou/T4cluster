#' CLUES algorithm
#' 
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param nbdsize size of the nearest neighborhood search for local shrinking (default: 100).
#' @param criterion type of quality measure, either \code{"CH"} or \code{"Sil"}.
#' @param kmax the upper bound for the number of clusters for search over (default: \eqn{n/5}).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{the maximum number of iterations (default: 10).}
#' \item{abstol}{stopping criterion (default: 1e-4).}
#' }
#' 
#' @return a named list of S3 class \code{T4cluster} containing 
#' \describe{
#' \item{cluster}{a length-\eqn{n} vector of class labels (from \eqn{1:k}).} 
#' \item{coords}{an \eqn{(n\times p)} matrix of shrunk data coordinates.}
#' \item{algorithm}{name of the algorithm.}
#' }
#' 
#' @examples 
#' \donttest{
#' # -------------------------------------------------------------
#' #            clustering with 'SMILEY' dataset
#' # -------------------------------------------------------------
#' ## PREPARE
#' smiley = T4cluster::genSMILEY(sd=0.25)
#' dat_true  = smiley$data
#' lab_true  = smiley$label
#' 
#' ## RUN CLUES WITH DIFFERENT NBDSIZES
#' run010 = clues(dat_true, kmax=5, nbdsize=10)
#' run050 = clues(dat_true, kmax=5, nbdsize=50)
#' run100 = clues(dat_true, kmax=5, nbdsize=100)
#' 
#' ## COMPARE VISUALIZATIONS OF SHRUNK DATA
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(run010$coords, pch=19, col=run010$cluster, main="nnbd=10")
#' plot(run050$coords, pch=19, col=run050$cluster, main="nnbd=50")
#' plot(run100$coords, pch=19, col=run100$cluster, main="nnbd=100")
#' par(opar)
#' }
#' 
#' @keywords internal
#' @noRd
clues <- function(data, nbdsize=10, criterion=c("Sil","CH"), kmax=round(nrow(data)/5), ...){
  ## PREPARE : EXPLICIT INPUTS
  mydata    = prec_input_matrix(data)
  mynbdsize = max(2, round(nbdsize))
  mycrit    = match.arg(criterion)
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  if ("maxiter"%in%pnames){
    myiter = max(1, round(pars$maxiter))
  } else {
    myiter = 100
  }
  if ("abstol"%in%pnames){
    myabstol = max(sqrt(.Machine$double.eps), as.double(pars$abstol))
  } else {
    myabstol = 1e-4
  }
  
  ## COMPUTE : SHRINKING
  shrink_data = clues_update(mydata, mynbdsize, myiter, myabstol)
  
  ## COMPUTE : FIND THE OPTIMAL 
  n     = base::nrow(mydata)
  k_min = 2
  k_max = min(round(kmax), n-1)
  k_seq = seq(from=2, to=k_max, by=1)
  k_num = length(k_seq)
  
  list_score   = rep(0, k_num)
  list_cluster = list()
  for (i in 1:k_num){
    # apply k-means 
    list_cluster[[i]] = T4cluster::kmeans(shrink_data, k=(i+1))$cluster
    if (all(mycrit=="CH")){
      list_score[i] = T4cluster::quality.CH(shrink_data, list_cluster[[i]])
    } else {
      list_score[i] = T4cluster::quality.sil(shrink_data, list_cluster[[i]])
    }
  }
  
  ## SELECT ONE AND RETURN
  max_id = which.max(list_score)
  max_cluster = list_cluster[[max_id]]
  
  ## SELECT, WRAP, AND RUN
  output = list()
  output$cluster   = max_cluster
  output$coords    = shrink_data
  output$algorithm ="clues"
  return(structure(output, class="T4cluster"))
}


# aux ---------------------------------------------------------------------
#' @keywords internal
#' @noRd
clues_update <- function(data, nnbd, maxiter, abstol){
  nbdfunc  = utils::getFromNamespace("hidden_knn","maotai")
  n = base::nrow(data)
  p = base::ncol(data)
  
  data_old = data
  data_new = array(0,c(n,p))
  for (it in 1:maxiter){
    # nearest neighbor information
    idx_info = nbdfunc(data, nnbd)$nn.idx
    
    # update each element
    for (i in 1:n){
      data_new[i,] = clue_median(data[as.vector(idx_info[i,]),])
    }
    
    # compute the stopping criterion
    dstop = max(abs(data_new-data_old))
    if (dstop < abstol){
      break
    } else {
      data_old <- data_new
    }
  }
  
  # return output
  return(data_old)
}
#' @keywords internal
#' @noRd
clue_median <- function(subdata){
  n = base::nrow(subdata)
  p = base::ncol(subdata)
  
  output = rep(0,p)
  for (i in 1:p){
    tgtvec = as.vector(subdata[,i])
    output[i] = stats::median(tgtvec)
  }
  return(output)
}