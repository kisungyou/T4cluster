## AUXILIARY FUNCTIONS
#  (01) prec_input_matrix : return output as row-stacked matrix
#       prec_input_dist   : return matrix of dist object
#       prec_input_sphere : return output as row-normalized matrix
#  (02) soc_preproc       : SOC-type algorithm preprocessing
#  (03) gmm_check_alldiag : check whether all covariances are diagonal or not
#  (04) gmm_name_segment  : extract the first 3 letters of algorithm name
#  (05) gmm_check_list    : check list of gmm objects
#  (06) gmm_barycenter    : given multiple models + support, compute opt.weight - obsolete
#  (07) prec_twolabel     : check the labels for 'mclustcomp' import

# (01) prec_input_matrix & prec_input_dist --------------------------------
#' @keywords internal
#' @noRd
prec_input_matrix <- function(x){
  if (is.vector(x)){
    output = matrix(x, ncol=1)
    return(output)
  } else {
    if (is.matrix(x)){
      return(x)
    } else {
      stop("* input should be either a vector or a matrix.")  
    }
  }
}
#' @keywords internal
#' @noRd
prec_input_sphere <- function(x){
  if (!inherits(x, "matrix")){
    stop("* input should be a matrix.")
  }
  return(x/sqrt(rowSums(x^2)))
}
#' @keywords internal
#' @noRd
prec_input_dist <- function(x){
  if (inherits(x, "dist")){
    return(as.matrix(x))
  } else if (is.vector(x)){
    x  = matrix(x, ncol=1)
    dx = cpp_pdist(x, 2)
    return(dx)
  } else if (is.matrix(x)){
    dx = cpp_pdist(x, 2)
    return(dx)
  } else {
    stop("* T4cluster : input should be either 'vector','matrix', or 'dist' object.")
  }
}



# (02) soc_preproc : SOC-type algorithm preprocessing ---------------------
#' @keywords internal
#' @noRd
soc_preproc <- function(partitions, fname){
  if (is.list(partitions)){
    M = length(partitions)
    if (length(unique(unlist(lapply(partitions, length)))) > 1){
      stop(paste0("* ",fname," : when a list is given, all element vectors should be of same length 
                  since they must be clustering vectors on the same-size data."))
    }
    N = length(partitions[[1]])
    clmat = array(0,c(M,N))
    for (m in 1:M){
      clmat[m,] = as.integer(as.factor(as.vector(partitions[[m]])))
    }
  } else {
    M = nrow(partitions)
    N = ncol(partitions)
    clmat = array(0,c(M,N))
    for (m in 1:M){
      clmat[m,] = as.integer(as.factor(as.vector(partitions[m,])))
    }
  }
  return(clmat)
}


# (03) gmm_check_alldiag --------------------------------------------------
#' @keywords internal
#' @noRd
gmm_check_alldiag <- function(covariance){
  k = dim(covariance)[3]
  for (i in 1:k){
    tgt = covariance[,,i]
    if (any(as.vector(tgt[upper.tri(tgt)])!=0)){
      return(FALSE)
    }
  }
  return(TRUE)
}

# (04) gmm_name_segment ---------------------------------------------------
#' @keywords internal
#' @noRd
gmm_name_segment <- function(strname){
  letter3 = paste(unlist(strsplit(strname,""))[1:3], collapse = "")
  return(letter3)
}


# (05) gmm_check_list : check list of gmm objects -------------------------
#' @keywords internal
#' @noRd
gmm_check_list <- function(gmmlist){
  if (!is.list(gmmlist)){
    return(FALSE)
  } else {
    K  = length(gmmlist)
    cc = rep(FALSE, K)
    for (k in 1:K){
      tgtgmm = gmmlist[[k]]
      cond1  = (inherits(tgtgmm, "T4cluster"))
      cond2  = ("algorithm"%in%names(tgtgmm))
      cond3  = (all(gmm_name_segment(tgtgmm$algorithm)=="gmm"))
      if (cond1&&cond2&&cond3){
        cc[k] = TRUE
      }
    }
    if (all(cc==TRUE)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

# (06) gmm_barycenter -----------------------------------------------------
#' #' @keywords internal
#' #' @noRd
#' gmm_barycenter <- function(gmmquery, gmmlist, lambda, maxiter, abstol){
#'   # preparation
#'   mean     = gmmquery$mean
#'   variance = gmmquery$variance
#'   
#'   N = length(gmmlist)
#'   K = base::nrow(mean)
#'   par_listdxy  <- list() # W_2 distances for sets of atoms
#'   par_marginal <- list() # marginal density for each gmm object
#'   par_weights  <- rep(1/N, N) # uniform weights across models
#'   par_iter  <- round(maxiter)    # Barycenter : iteration
#'   par_tol   <- as.double(abstol) # Barycenter : absolute tolerance
#'   par_lbd   <- as.double(lambda) # Barycenter : regularization parameters
#'   par_p     <- 2                 # Barycenter : we are working with W2
#'   # par_init  <- rep(1/K, K)       # Barycenter : initial weight
#'   par_init  <- as.vector(gmmquery$weight)
#'   printer   <- FALSE
#'   
#'   par_iter <- max(200, par_iter)
#'   par_tol  <- max(1e-10, par_tol)
#'   par_lbd  <- max(1e-10, par_lbd)
#'   
#'   # assignment
#'   for (n in 1:N){
#'     ngmm = gmmlist[[n]]
#'     par_marginal[[n]] = as.vector(ngmm$weight)
#'     par_listdxy[[n]]  = cpp_gmmdist_base(mean, variance, ngmm$mean, ngmm$variance, "wass2")
#'   }
#'   
#'   # main computation
#'   output = cpp_barybregman15(par_listdxy, par_marginal, par_weights,
#'                              par_p, par_lbd, par_iter, par_tol, printer, par_init)
#'   return(as.vector(output))
#' }

# (07) prec_twolabel     : check the labels for 'mclustcomp' import -------
#' @keywords internal
#' @noRd
prec_twolabel <- function(x, y, fname){
  if ((!is.vector(x))||(!is.vector(y))){
    stop(paste0("* ",fname," : input 'x' and 'y' should both be a vector of class labels."))
  }
  if (length(x)!=length(y)){
    stop(paste0("* ",fname," : two vectors should be of same size."))
  }
}


