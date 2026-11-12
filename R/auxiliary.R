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
#  (08) check_clist       : a list of consensus matrices

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
prec_input_sphere <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x)) {
    stop("data must be a real numeric matrix.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 2L) {
    stop("data must have at least one row and two columns.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("data must contain only finite values.", call. = FALSE)
  }

  # Scale before squaring: finite nonzero rows can be very large or very small.
  row_scale <- apply(abs(x), 1L, max)
  zero_rows <- which(row_scale == 0)
  if (length(zero_rows)) {
    stop(
      paste0("Zero rows have no spherical direction. Row indices include: ",
             paste(utils::head(zero_rows, 8L), collapse = ", "), "."),
      call. = FALSE
    )
  }
  scaled <- sweep(x, 1L, row_scale, "/")
  output <- sweep(scaled, 1L, sqrt(rowSums(scaled * scaled)), "/")
  storage.mode(output) <- "double"
  output
}

#' @keywords internal
#' @noRd
prec_sphere_controls <- function(k, n, dots, intrinsic = FALSE) {
  scalar_integer <- function(x, name, upper = .Machine$integer.max) {
    valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L &&
      is.finite(x) && x >= 1 && x <= upper && x == floor(x)
    if (!valid) {
      stop(sprintf("%s must be one finite integer between 1 and %s.",
                   name, format(upper, scientific = FALSE)), call. = FALSE)
    }
    as.integer(x)
  }
  scalar_positive <- function(x, name) {
    valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L &&
      is.finite(x) && x > 0
    if (!valid) {
      stop(paste0(name, " must be one finite positive number."), call. = FALSE)
    }
    as.double(x)
  }

  allowed <- c("init", "maxiter", "abstol", "verbose")
  if (intrinsic) allowed <- c(allowed, "mean.maxiter", "mean.tol")
  dot_names <- names(dots)
  if (length(dots)) {
    if (is.null(dot_names) || anyNA(dot_names) || any(!nzchar(dot_names))) {
      stop("All arguments in ... must be named.", call. = FALSE)
    }
    if (anyDuplicated(dot_names)) {
      stop("Arguments in ... must not be duplicated.", call. = FALSE)
    }
    unknown <- setdiff(dot_names, allowed)
    if (length(unknown)) {
      stop(paste0("Unknown argument(s): ", paste(unknown, collapse = ", "),
                  ". Allowed arguments: ", paste(allowed, collapse = ", "), "."),
           call. = FALSE)
    }
  }
  get_option <- function(name, default) {
    if (name %in% dot_names) dots[[name]] else default
  }
  init <- get_option("init", "kmeans")
  if (!is.character(init) || length(init) != 1L || is.na(init)) {
    stop("init must be one character value: 'kmeans' or 'gmm'.", call. = FALSE)
  }
  init <- match.arg(tolower(init), c("kmeans", "gmm"))
  verbose <- get_option("verbose", FALSE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  control <- list(
    k = scalar_integer(k, "k", upper = min(n, .Machine$integer.max)),
    init = init,
    maxiter = scalar_integer(get_option("maxiter", 10L), "maxiter"),
    abstol = scalar_positive(
      get_option("abstol", sqrt(.Machine$double.eps)), "abstol"),
    verbose = verbose
  )
  if (intrinsic) {
    control$mean.maxiter <- scalar_integer(
      get_option("mean.maxiter", 100L), "mean.maxiter")
    control$mean.tol <- scalar_positive(get_option("mean.tol", 1e-8), "mean.tol")
  }
  control
}

#' @keywords internal
#' @noRd
prec_sphere_result <- function(cpprun, algorithm, control) {
  output <- list(
    cluster = as.integer(cpprun$cluster) + 1L,
    cost = as.double(cpprun$cost),
    means = cpprun$means,
    algorithm = algorithm
  )
  reserved <- c("cluster", "cost", "means", "algorithm", "control")
  extra_names <- setdiff(names(cpprun), reserved)
  if (length(extra_names)) output[extra_names] <- cpprun[extra_names]
  output$control <- control
  structure(output, class = "T4cluster")
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


# (08) check_clist --------------------------------------------------------
#' @keywords internal
#' @noRd
check_clist <- function(clist, fname){
  if (!is.list(clist)){
    stop(paste0("* ",fname," : 'clist' must be a list."))
  }
  vec_class = unlist(base::lapply(clist, inherits, "T4cluster:consensus"))
  if (any(!vec_class)){
    stop(paste0("* ",fname," : 'clist' must consist of 'T4cluster:consensus' objects."))
  }
  if (length(unique(unlist(base::lapply(clist, nrow))))!=1){
    stop(paste0("* ",fname," : all consensus matrices should be of same size."))
  }
}
