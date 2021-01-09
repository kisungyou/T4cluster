#' Generate Line and Plane Example with Fixed Number of Components
#' 
#' This function generates a toy example of 'line and plane' data in \eqn{R^3} that 
#' data are generated from a mixture of lines (one-dimensional) planes (two-dimensional).
#' The number of line- and plane-components are explicitly set by the user for flexible testing.
#' 
#' @param n the number of data points for each line and plane.
#' @param nl the number of line components.
#' @param np the number of plane components.
#' @param iso.var degree of isotropic variance.
#' 
#' @return a named list containing with \eqn{m = n\times(nl+np)}:\describe{
#' \item{data}{an \eqn{(m\times 3)} data matrix.}
#' \item{class}{length-\eqn{m} vector for class label.}
#' \item{dimension}{length-\eqn{m} vector of corresponding dimension from which an observation is created.}
#' }
#' 
#' @examples 
#' ## test for visualization
#' set.seed(10)
#' tester = genLP(n=100, nl=1, np=2, iso.var=0.1)
#' data   = tester$data
#' label  = tester$class
#' 
#' ## do PCA for data reduction
#' proj = base::eigen(stats::cov(data))$vectors[,1:2]
#' dat2 = data%*%proj
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.5,col=label,main="PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.5,col=label,main="Axis 1 vs 2")
#' plot(data[,1],data[,3],pch=19,cex=0.5,col=label,main="Axis 1 vs 3")
#' plot(data[,2],data[,3],pch=19,cex=0.5,col=label,main="Axis 2 vs 3")
#' par(opar)
#' 
#' \dontrun{
#' ## visualize in 3d
#' x11()
#' scatterplot3d::scatterplot3d(x=data, pch=19, cex.symbols=0.5, color=label)
#' }
#' 
#' @concept utility
#' @export
genLP <- function(n=100, nl=1, np=1, iso.var=0.1){
  # parameters
  m=3
  K=round(nl)+round(np)
  n=round(n)
  isotropic.var=as.double(iso.var)
  k.true= c(rep(1, round(nl)), rep(2, round(np)))
  
  #J.true=2
  #Generate the subspaces
  Utrue.list = list()
  for(k in 1:K){
    Utrue.list[[k]] = rstiefel::rustiefel(m=m,R=k.true[k])
  }
  NU.true = list()
  for(k in 1:K){
    NU.true[[k]]=MASS::Null(Utrue.list[[k]])
  }
  PNU.list = list()
  for(k in 1:K){
    PNU.list[[k]] = NU.true[[k]]%*%t(NU.true[[k]])
  }
  #generate means in subspace coordinates
  mutrue.list = list()
  for(k in 1:K){
    mutrue.list[[k]] = rnorm(k.true[k])
  }
  
  #generate the residual space noise level
  phitrue.list = rep(10,K)
  sigmatrue.list = rep(isotropic.var,K)
  #generate the subspace variances
  sigma0.true.list = list()
  for(k in 1:K){
    sigma0.true.list[[k]]=runif(k.true[k],isotropic.var,5.1)
  }
  Sigma0.true.list = list()
  for(k in 1:K){
    Sigma0.true.list[[k]]=diag(sigma0.true.list[[k]],k.true[k])
  }
  
  
  #generate the euclidean space coordinate mean vector theta
  theta.true.list = list()
  for(k in 1:K){
    theta.true.list[[k]] = PNU.list[[k]]%*%rnorm(m)
  }
  X = c()
  label = c()
  dimension = c()
  for(k in 1:K){
    X = rbind(X,MASS::mvrnorm(n,mu=Utrue.list[[k]]%*%mutrue.list[[k]]+theta.true.list[[k]],
                              Sigma = sigmatrue.list[k]^2*diag(m)/10+Utrue.list[[k]]%*%(Sigma0.true.list[[k]]-sigmatrue.list[k]*diag(k.true[k]))%*%t(Utrue.list[[k]])))
    label = c(label, rep(k,n))
    dimension = c(dimension, rep(k.true[k], n))
  }
  # scatterplot3d(x=X[,1], y=X[,2], z=X[,3], pch = 19, color = label)
  
  
  # return
  output = list()
  output$data  = X
  output$class = label
  output$dimension = dimension
  return(output)
}