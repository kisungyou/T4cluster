#' Generate from Three 5-dimensional Subspaces in 200-dimensional space.
#' 
#' @param n the number of data points sampled from each subspace (default: 50).
#' @param var degree of Gaussian noise (default: 0.3).
#' 
#' @return a named list containing with :\describe{
#' \item{data}{an \eqn{(3*n\times 3)} data matrix.}
#' \item{class}{length-\eqn{3*n} vector for class label.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## a toy example
#' tester = gen3S(n=100)
#' data   = tester$data
#' label  = tester$class
#' }
#' 
#' @references 
#' \insertRef{wang_efficient_2011}{T4cluster}
#' 
#' @concept utility
#' @export
gen3S <- function(n=50, var=0.3){
  ## PRELIMINARY
  par_n  = max(2, round(n))
  par_sd = base::sqrt(max(10*.Machine$double.eps, as.double(var)))
  
  ## LET'S GENERATE !
  rot1 = qr.Q(qr(matrix(stats::rnorm(200*200),ncol=200)))
  rot2 = qr.Q(qr(matrix(stats::rnorm(200*200),ncol=200)))
  
  U1 = qr.Q(qr(matrix(stats::rnorm(200*5),ncol=5))) 
  U2 = rot1%*%U1
  U3 = rot2%*%U2
  
  Q1 = matrix(stats::runif(5*par_n), ncol = par_n)
  Q2 = matrix(stats::runif(5*par_n), ncol = par_n)
  Q3 = matrix(stats::runif(5*par_n), ncol = par_n)
  
  X1 = U1%*%Q1 + matrix(stats::rnorm(200*par_n, sd=par_sd), ncol=par_n)
  X2 = U2%*%Q2 + matrix(stats::rnorm(200*par_n, sd=par_sd), ncol=par_n)
  X3 = U3%*%Q3 + matrix(stats::rnorm(200*par_n, sd=par_sd), ncol=par_n)
  
  # return
  output = list()
  output$data  = t(cbind(X1, X2, X3))
  output$class = rep(1:3, each=par_n)
  return(output)
}