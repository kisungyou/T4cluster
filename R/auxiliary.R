## AUXILIARY FUNCTIONS
#  (01) prec_input_matrix : return output as row-stacked matrix
#       prec_input_dist   : return matrix of dist object
#  (02) soc_preproc       : SOC-type algorithm preprocessing


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

#  (02) soc_preproc       : SOC-type algorithm preprocessing -------------------
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