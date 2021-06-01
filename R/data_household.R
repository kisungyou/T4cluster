#' Load 'household' data
#' 
#' The data is taken from \pkg{HSAUR3} package's \code{household} data. We use 
#' housing, service, and food variables and normalize them to be unit-norm so 
#' that each observation is projected onto the 2-dimensional sphere. The data 
#' consists of 20 males and 20 females and has been used for clustering 
#' on the unit hypersphere.
#' 
#' @usage data(household)
#' 
#' @examples
#' \donttest{
#' ## Load the data
#' data(household, package="T4cluster")
#' 
#' ## Visualize the data in pairs
#' opar <- par(no.readonly=TRUE)
#' scatterplot3d::scatterplot3d(household$data, color=rep(c("red","blue"), each=20), 
#'               pch=19, main="household expenditure on the 2-dimensional sphere",
#'               xlim=c(0,1.2), ylim=c(0,1.2), zlim=c(0,1.2), angle=45)
#' par(opar)
#' }
#' 
#' @format a named list containing\describe{
#' \item{data}{an \eqn{(n\times 3)} data matrix whose rows are unit-norm.}
#' \item{gender}{a length-\eqn{n} factor for class label.}
#' }
#' 
#' @seealso \code{\link[HSAUR3]{household}}
#' @concept data
"household"