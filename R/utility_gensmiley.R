#' Generate SMILEY Data
#' 
#' Creates a smiley-face data in \eqn{\mathbf{R}^2}. This function is a modification 
#' of \pkg{mlbench}'s \code{mlbench.smiley} function.
#' 
#' @param n number of samples to be generated.
#' @param sd additive Gaussian noise level.
#' 
#' @return a list containing\describe{
#' \item{data}{an \eqn{(n\times 2)} data matrix.}
#' \item{label}{a length-\eqn{n} vector(factor) for class labels.}
#' } 
#' 
#' @examples 
#' ## Generate SMILEY Data with Difference Noise Levels
#' s10 = gensmiley(200, sd=0.1)
#' s25 = gensmiley(200, sd=0.25)
#' s50 = gensmiley(200, sd=0.5)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(s10$data, col=s10$label, pch=19, main="sd=0.10")
#' plot(s25$data, col=s25$label, pch=19, main="sd=0.25")
#' plot(s50$data, col=s50$label, pch=19, main="sd=0.50")
#' par(opar)
#' 
#' @seealso \code{\link[mlbench]{mlbench.smiley}}
#' @concept utility
#' @export
gensmiley <- function(n=496, sd=0.1){
  # parameters
  myn  <- max(50, round(n))
  mysd <- as.double(sd)
  
  # sizes of each object
  n_eye   <- round(n/6)
  n_nose  <- round(n/4)
  n_mouth <- (myn-(2*n_eye)-n_nose)
  
  # data generation
  dat_eye1 <- cbind(rnorm(n_eye, -.8, sd), rnorm(n_eye, 1, sd))
  dat_eye2 <- cbind(rnorm(n_eye,  .8, sd), rnorm(n_eye, 1, sd))
  
  dat_nose <- cbind(runif(n_nose, -.2, .2), runif(n_nose, 0, .75))
  dat_nose[,1] <- dat_nose[,1]*(1-dat_nose[,2])
  
  dat_mouth <- sort(runif(n_mouth, -1, 1))
  dat_mouth <- cbind(dat_mouth, (dat_mouth^2) - 1 + rnorm(n_mouth, 0, sd))
  
  output = list()
  output$data  = as.matrix(rbind(dat_eye1, dat_eye2, dat_nose, dat_mouth))
  output$label = as.factor(c(rep(1,n_eye), rep(2,n_eye), rep(3,n_nose), rep(4,n_mouth)))
  return(output)
}