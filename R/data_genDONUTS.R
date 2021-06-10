#' Generate Nested Donuts
#' 
#' It generates nested \emph{donuts}, which are just hollow circles. For flexible 
#' testing, the parameter \code{k} controls the number of circles of varying 
#' radii where \code{n} controls the number of observations for each circle. 
#' 
#' @param n the number of data points for each hollow circle (default: 50).
#' @param k the number of circles (default: 2).
#' @param sd magnitude of white noise (default: 0.1).
#' 
#' @return a named list containing with \eqn{m = nk}:\describe{
#' \item{data}{an \eqn{(m\times 2)} data matrix.}
#' \item{label}{a length-\eqn{m} vector(factor) for class labels.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## generate data
#' donut2 = genDONUTS(k=2)
#' donut3 = genDONUTS(k=3)
#' donut4 = genDONUTS(k=4)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(donut2$data, col=donut2$label, pch=19, main="k=2")
#' plot(donut3$data, col=donut3$label, pch=19, main="k=3")
#' plot(donut4$data, col=donut4$label, pch=19, main="k=4")
#' par(opar)
#' }
#' 
#' @concept data
#' @export
genDONUTS <- function(n=50, k=2, sd=0.1){
  # parameters
  myn  = max(1, round(n))
  myk  = max(2, round(k))
  mysd = max(.Machine$double.eps, as.double(sd))
  
  # generate data
  dat = c()
  for (i in 1:myk){
    theta.i = stats::runif(myn, min=0, max=2*pi)
    data.i  = cbind(i*cos(theta.i), i*sin(theta.i)) + matrix(stats::rnorm(50*2, sd=mysd),ncol=2)
    dat = rbind(dat, data.i)
  }
  lab = rep(1:myk, each=myn)
  

  # return
  output = list()
  output$data  = dat
  output$label = factor(lab)
  return(output)
}