## AUXILIARY FUNCTIONS
#  (01) prec_input_matrix : return output as row-stacked matrix
#       prec_input_dist   : return matrix of dist object


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
