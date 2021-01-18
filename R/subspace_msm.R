#' Bayesian Mixture of Subspaces of Different Dimensions 
#' 
#' \code{msm} is a Bayesian model inferring mixtures of subspaces that are of possibly different dimensions. 
#' For simplicity, this function returns only a handful of information that are most important in 
#' representing the mixture model, including projection, location, and hard assignment parameters. 
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @param k the number of mixtures.
#' @param ... extra parameters including \describe{
#' \item{temperature}{temperature value for Gibbs posterior (default: 1e-6).}
#' \item{prop.var}{proposal variance parameter (default: 1.0).}
#' \item{iter}{the number of MCMC runs (default: 496).}
#' \item{burn.in}{burn-in for MCMC runs (default: iter/2).}
#' \item{thin}{interval for recording MCMC runs (default: 10).}
#' \item{print.progress}{a logical; \code{TRUE} to show completion of iterations by 10, \code{FALSE} otherwise (default: \code{FALSE}).}
#' }
#' 
#' @return a list whose elements are S3 class \code{"msm"} instances, which are also lists of following elements: \describe{
#' \item{P}{length-\code{k} list of projection matrices.}
#' \item{U}{length-\code{k} list of orthonormal basis.}
#' \item{theta}{length-\code{k} list of center locations of each mixture.}
#' \item{cluster}{length-\code{n} vector of cluster label.}
#' }
#' 
#' @examples
#' \donttest{
#' ## generate a toy example
#' set.seed(10)
#' tester = genLP(n=100, nl=2, np=1, iso.var=0.1)
#' data   = tester$data
#' label  = tester$class
#' 
#' ## do PCA for data reduction
#' proj = base::eigen(stats::cov(data))$vectors[,1:2]
#' dat2 = data%*%proj
#' 
#' ## run MSM algorithm with k=2, 3, and 4
#' maxiter = 500
#' output2 = msm(data, k=2, iter=maxiter)
#' output3 = msm(data, k=3, iter=maxiter)
#' output4 = msm(data, k=4, iter=maxiter)
#' 
#' ## extract final clustering information
#' nrec  = length(output2)
#' finc2 = output2[[nrec]]$cluster
#' finc3 = output3[[nrec]]$cluster
#' finc4 = output4[[nrec]]$cluster
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(3,4))
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc2+1,main="K=2:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc3+1,main="K=3:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc4+1,main="K=4:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(2,3)")
#' par(opar)
#' }
#' 
#' @concept subspace
#' @export
msm <- function(data, k=2, ...){
  ## PREPARE : EXPLICIT INPUTS
  X = prec_input_matrix(data)
  m = ncol(X)
  n = nrow(X)
  K = max(1, round(k))
  
  ## PREPARE : IMPLICIT INPUTS
  pars     = list(...)
  pnames   = names(pars)
  
  temperature = ifelse("temperature"%in%pnames, pars$temperature, 1e-6)
  talk = ifelse("print.progress"%in%pnames, pars$print.progress, FALSE)
  report = 10
  sub.met    = ifelse("prop.var"%in%pnames, pars$prop.var, 1.0)
  sub.met.sd = rep(sub.met, K)
  R = base::sample(1:(m-1), K, replace=TRUE)
  Urandom  = FALSE # local pca
  my.kappa = 10
  iter    = ifelse("iter"%in%pnames, pars$iter, 496)
  burn.in = ifelse("burn.in"%in%pnames, pars$burn.in, round(iter/2))
  thin    = ifelse("thin"%in%pnames, pars$thin, 10)
  
  ############# copied code
  #Intialize the list of K subspaces and projection matrices
  kmeans.init = stats::kmeans(X,K)
  U.list = list()
  P.list = list()
  if (Urandom){
    for(k in 1:K){
      unow        = rustiefel(m=m, R = R[k])
      U.list[[k]] = unow
      P.list[[k]] = (unow%*%t(unow))
    }  
  } else {
    for (k in 1:K){
      idnow = which(kmeans.init$cluster==k)
      U.list[[k]] = base::eigen(stats::cov(X[idnow,]))$vectors[,(1:R[k])]
      P.list[[k]] = kisung.outer(U.list[[k]])
    }
  }
  #Initialize theta list
  theta.list = list()
  for(k in 1:K){
    theta.list[[k]] =Null(U.list[[k]])%*%t(Null(U.list[[k]]))%*% kmeans.init$centers[k,]
  } 
  
  #Initialize theta storage
  theta.mat = vector("list",iter)
  
  #Intialize subspace storage. Each subspace iteration will store U.mat[[iter]] = U.list
  #Allowing each subspace to be accessed as U.mat[[iter]][[k]]
  U.mat = vector("list",iter)
  P.mat = vector("list",iter)
  
  #Initialize the distance matrix to be stored every iteration
  distance = array(0,c(K,n)) 
  
  #Initialize storage for the latent normal and sphere walks 
  z.store = array(0,dim = c(iter,K,m*(m+1)/2))
  s.store = array(0,dim = c(iter,K,m*(m+1)/2))
  
  #initialize a random normal location
  z = matrix(0,K,m*(m+1)/2)
  s = matrix(0,K,m*(m+1)/2)
  for (k in 1:K){
    temp = conway.sphere.step(z=rep(0,m*(m+1)/2),sd=1,m=m)
    z[k,]=temp$z
    s[k,]=temp$s
  }
  
  #Initialize acceptance counter
  accept = rep(0,K)
  
  #Get loss for initial estimates
  curr.lossclus = fast.log.loss(x=X, P.list = P.list,mu = theta.list,temperature=temperature)
  curr.loss = curr.lossclus$loss
  curr.clus = curr.lossclus$clus
  
  #initialize and store the clustering
  lat.clus = which(rmultinom(n=n, size = 1, prob=rep(1/K,K))==1,arr.ind=T)[,1]
  lat.clus.mat = matrix(0,nrow=iter,ncol=n)
  pi.mat = matrix(0,n,K)
  
  #Intialize the multinomial weights
  pi.mat = matrix(1/K,nrow=n,ncol=K)
  dir.prior = rep(1/K,K)
  n.vec = rep(0,K)
  r0 = rep(1,K)
  
  #set up tuning parameter
  tune = 0 
  tune.accept = rep(0,K)
  record.clus = list()
  for(i in 1:iter){
    # print(paste('iter is ',i))
    tune = tune+1
    if(talk){
      if(i%%report==0){
        print(paste('* msm : iteration ', i,'/',iter,' complete at ', Sys.time(),sep=""))
        
      }
    }
    
    
    #For each component, generate a proposal subspace, and then
    #accept or reject it based on the gibbs posterior
    for(k in 1:K){
      #Get proposal projection matrix
      #print('get proposal')
      proposal = conway.step(z=z[k,],sd=sub.met.sd[k],m=m)
      #print('get Subspace')
      prop.sub = con2sub(P=unembed(proposal$s),return.proj = F)
      #restrict samples to m/2 ([m+1]/2 if m is odd ) to stay on lower half of sphere
      #print('Restrict')
      if(is.matrix(prop.sub)){
        if(dim(prop.sub)[2]>ceiling(m/2)){
          if(dim(prop.sub)[2]==m){
            #print('Proposal full dimension, replace with 1 dimension')
            prop.sub= rustiefel(R=1,m=m)
            proposal$z=embed(prop.sub,subspace=T)
          }else{
            #print('Proposal not full dimension')
            prop.sub = Null(prop.sub)
            proposal$z = embed(prop.sub,subspace=T)
          }
        }
      }
      # print('Get projection')
      prop.proj = prop.sub%*%t(prop.sub)
      
      #Set the proposal list
      prop.P.list = P.list
      prop.P.list[[k]] = prop.proj
      
      #Choose an appropriate theta in the null space
      prop.null = Null(prop.sub)
      prop.nullproj = prop.null%*%t(prop.null)
      prop.theta.list = theta.list
      prop.theta.list[[k]] = prop.nullproj%*%theta.list[[k]]
      for(l in 1:K){
        if(is.null(dim(prop.theta.list[[l]]))){
          prop.theta.list[[l]]=matrix(prop.theta.list[[l]],nrow=m,ncol=1)
        }
      }
      
      #Get the loss of the proposal
      prop.lossclus = fast.log.loss(x=X,P.list = prop.P.list, mu = prop.theta.list,temperature=temperature)
      prop.loss = prop.lossclus$loss
      #The coin toss
      toss = log(runif(1,0,1))
      diff.loss = prop.loss - curr.loss
      if(toss<diff.loss){
        accept[k] = accept[k] +1
        tune.accept[k] = tune.accept[k]+1
        #int.acc = int.acc +1
        P.list[[k]] = prop.proj
        U.list[[k]] = (prop.sub) #--------------------------------
        theta.list[[k]] = prop.theta.list[[k]]
        z[k,] =proposal$z 
        curr.loss = prop.loss
        curr.clus = prop.lossclus$clus
      }
    }
    
    # ## ADD : I want no cluster to be empty
    if (length(unique(curr.clus))<K){
      add.lossclus = fast.log.loss(x=X,P.list=P.list,mu=theta.list,temperature=temperature)
      add.wi       = mygibbs.step.b(t(add.lossclus$dist), kappa=my.kappa)
      curr.clus    = mygibbs.step.c(add.wi)  
    }
    record.clus[[i]] = curr.clus
    
    #Set new means based on closest clustering
    for(k in 1:K){
      idnow = which(curr.clus==k)
      if (length(idnow)==1){
        theta.temp = X[idnow,]
      } else if (length(idnow)>1){
        theta.temp = apply(X[curr.clus==k,],2,mean);  
      } else {
        stop("* there is no cluster.")
      }
      theta.list[[k]]=theta.temp
      theta.list[[k]] = Null(U.list[[k]])%*%t(Null(U.list[[k]]))%*%theta.list[[k]]  
    }
    
    #increase or decrease variance to adjust acceptance rates
    if(tune == 100){
      for(k in 1:K){
        if(tune.accept[k]<10){
          sub.met.sd[k] = .1*sub.met.sd[k]
        }else{
          if(tune.accept[k]<30){
            sub.met.sd[k]=.5*sub.met.sd[k]
          }else{
            if(tune.accept[k]<60){
              sub.met.sd[k]=2*sub.met.sd[k]
            }else{
              if(tune.accept[k]<90){
                sub.met.sd[k] = 5*sub.met.sd[k]
              }else{
                sub.met.sd[k]=10*sub.met.sd[k]
              }
            }
          }
        }
        tune.accept[k]=0
      }
      tune=0
    }
    
    
    #For storage at the end
    P.mat[[i]] = P.list
    U.mat[[i]] = U.list
    theta.mat[[i]]=theta.list
  }
  
  # return output
  output = list()
  for (i in 1:iter){
    iterate = list()
    iterate$P = P.mat[[i]]
    iterate$U = U.mat[[i]]
    iterate$theta = theta.mat[[i]]
    iterate$cluster = record.clus[[i]]
    output[[i]] = structure(iterate, class="msm")
  }
  recording = seq(from=round(burn.in+1), to=iter, by = round(thin))
  return(output[recording])
}


# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
kisung.outer <- function(A){
  return(A%*%t(A))
}
#' @keywords internal
#' @noRd
Trace = function(X){
  return(sum(diag(X)))
}
#' @keywords internal
#' @noRd
embed = function(P, subspace = F){
  #returns the conway embedding of the projection matrix P if subspace =F, 
  #or the subspace P if subspace =T
  m=dim(P)[1]
  if(subspace == T){
    P = P%*%t(P)
  }
  p=c()
  for(i in 1:m){
    p=c(p,P[i:m,i])
  }  
  return(p)
}
#' @keywords internal
distance = function(U,V,subspace=T){
  #Find distance between m-dimensional subspaces
  m=dim(U)[1];
  d = sqrt(m)
  if(subspace==T){
    P1 = U%*%t(U);
    P2 = V%*%t(V);
  }else{
    P1=U;
    P2=V;
  }
  PZ = diag(.5,m)
  p1 = embed(P1)
  p2 = embed(P2)
  pz = embed(PZ)
  angle = acos(sum((p1-pz)*(p2-pz))/(sqrt(sum((p1-pz)^2))*sqrt(sum((p2-pz)^2))))
  return(d*angle/2)
}
#' @keywords internal
#' @noRd
con.distance = function(U,V,subspace=T){
  #Find distance between ambient m-dimensional subspaces
  #using the conway distance
  #aka projection distance
  m=dim(U)[1];
  d = sqrt(m)
  if(subspace==F){
    temp1 = eigen(U);
    U = temp1$vectors[,temp1$values>10^(-12)]
    temp2 = eigen(V);
    V = temp2$vectors[,temp2$values>10^(-12)]
  }
  sigma = svd(t(U)%*%V)$d
  sigma[sigma>1]=1
  princ = acos(sigma)
  return(sqrt(sum((sin(princ))^2)))
  
}
#' @keywords internal
#' @noRd
unembed = function(p){
  #Given an embedding of a subspace p,
  #return the Projection matrix P
  m = (sqrt(8*length(p)+1)-1)/2
  P = matrix(0,nrow=m,ncol=m);
  take = m:1
  place = 1
  for(i in 1:m){
    P[i:m,i]=p[place:(place+take[i]-1)]
    place=place+take[i]
  }
  P = P + t(P)-diag(diag(P),m)
  return(P)
}
#' @keywords internal
#' @noRd
con2sub = function(P,d=Inf,return.proj = T){
  #go from the embededd space to either the 
  #projection matrix or the subspace depending
  #on return.proj=T
  
  #get dimensions
  m = dim(P)[1]
  
  if(d == Inf){
    
    prj.rank = sum(diag(P))
    
    if(prj.rank%%1>.5){
      d=ceiling(prj.rank)
    }else{
      d=floor(prj.rank)
    }
    
  }
  temp = eigen(P)
  U = temp$vectors[,1:d]
  if(return.proj ==T){
    return(U%*%t(U))
  }else{
    return(U)
  }
}
#' @keywords internal
#' @noRd
runif.sphere = function(p=3){
  #sample a p-dimensional vector from the unit sphere in R^p
  v = rnorm(p,0,1)
  v=v/sqrt(sum(v^2))
  return(v)
}
#' @keywords internal
#' @noRd
runif.conway = function(n=1,m=3){
  p  = m*(m+1)/2
  r  = sqrt(m)/2
  PZ = diag(rep(.5,m),m)
  pz = embed(PZ)
  v  = array(0,c(n,p))
  for (i in 1:n){
    v[i,] = as.vector(r*runif.sphere(p)+pz)
  }
  if(n==1){
    return(as.numeric(v))
  }else{
    return(v)
  }
}
#' @keywords internal
#' @noRd
conway.spherize = function(z,m=3){
  #Renormalize the point Z onto the conway sphere
  #by first mapping it to the unit sphere, and then 
  #translating and scaling to the Conway sphere
  z = z/sqrt(sum(z^2))
  pz = embed(diag(rep(.5,m),m))
  r = sqrt(m)/2
  return(z*r+pz)
}
#' @keywords internal
#' @noRd
embed.rustiefel = function(n,m,R){
  X = matrix(0,nrow=n,ncol=m*(m+1)/2)
  for(i in 1:n){
    X[i,]= embed(rustiefel(m=m,R=R),subspace=T)
  }
  return(X)
}
#' @keywords internal
#' @noRd
normal.walk = function(n=1000,m=3){
  v = matrix(0,nrow=n,ncol=m)
  v[1,]=rnorm(m)
  for( i in 2:n){
    v[i,]=rnorm(m)+v[i-1,]
  }
  return(v)
}
#' @keywords internal
#' @noRd
sphere.walk = function(n=1000,m=3){
  #generates a random walk on the unit sphere of length n
  v = matrix(0,nrow=n,ncol=m)
  s = matrix(0,nrow=n,ncol=m)
  v[1,]=rnorm(m)
  s[1,]=v[1,]/sum(v[1,]^2)
  for(i in 2:n){
    v[i,]=rnorm(m)+v[i-1,]
    s[i,] = v[i,]/sum(v[i,]^2)
  }
  return(list('v'=v,'s'=s))
}
#' @keywords internal
sphere.step = function(z,sd=1,m=3){
  #given a specific location on the latent walk,
  #talk a step and spherize it. Return the new
  #location
  p = m*(m+1)/2
  v = rnorm(p,mean=z,sd=sd)
  s = v/sqrt(sum(v^2))
  return(list('z'=v,'s'=s))
}
#' @keywords internal
conway.sphere.step = function(z, sd=1,m=3){
  #Given a specific location of the latent walk
  #take a step and spherize it. Return both the new location
  #of the latent walk and the spherized step. Now the 
  #'s' is given on the Conway sphere in m(m+1)/2 space
  r= sqrt(m)/2
  PZ = diag(rep(.5,m),m)
  pz = embed(PZ)
  step = sphere.step(z,sd=sd,m=m)
  s=r*step$s+pz
  return(list('s'=s,'z'=step$z))
}
#' @keywords internal
conway.step = function(z,sd=1,m=3){
  #Given a location on the Conway sphere z,
  #a standard deviation sd, and the ambient dimension m
  # take a random step on the sphere, and return the sphere location
  r = sqrt(m)/2;
  p = m*(m+1)/2
  PZ = diag(rep(.5,m),m);
  pz= embed(PZ)
  tangent.step = rnorm(p,z,sd)
  tangent.step = tangent.step/sqrt(sum(tangent.step^2))
  step = tangent.step*r+pz
  return(list('s'=step,'z'=tangent.step))
}
#' @keywords internal
gibbs.loss.prj=function(x , P, mu = NULL, subspace = F){
  #P is one a  m x m projection matrices (or a subspace)
  # x is the set of n observations, so is n x m. Finds the distance
  #between theclosest point on the subspace P and x
  n = dim(x)[1]
  if(subspace==T){
    # P=P%*%P^T
    P = P%*%t(P)
  }
  m = dim(P)[1]
  d = sum(diag(P))
  
  if(is.null(mu)){
    mu = rep(0,m)
  } else {
    mu = as.vector(mu)
  }
  
  #Find the distance between the projection of the point and the point itself for each 
  #observation
  PIm  = P-diag(m)
  dist = rep(0, n)
  for(i in 1:n){
    xidiff  = x[i,]-mu
    dist[i] = d*sqrt(sum(m^2*(as.vector(PIm%*%xidiff))^2))
  }
  
  return(dist)
  
}
#' @keywords internal
log.gibbs.loss = function(x, P.list, mu = Inf,temperature=1, subspace =F){
  
  #P.list is a list of projection matrices, indexed by P.list[[k]]
  #mu is the list of theta, where each can be accessed as mu[[k]]
  #given a set of observations x
  K = length(P.list)
  n = dim(x)[1]
  dist = matrix(0,nrow = K, ncol = n)
  for(k in 1:K){
    #print('log gibbs loss')
    # print(k)
    dist[k,] = gibbs.loss.prj(x, P.list[[k]],mu = mu[[k]],subspace=subspace)
  }
  mindist = apply(dist,2,min)
  whichmindist = apply(dist,2,which.min)
  loss = -1*n*temperature *sum(mindist)
  return(list('loss'=loss,'clus'=whichmindist))
}

#' @keywords internal
#' @noRd
fast.log.loss = function(x, P.list, mu, temperature = 1, subspace =F){
  #P.list is a list of projection matrices, indexed by P.list[[k]]
  #mu is the list of theta, where each can be accessed as mu[[k]]
  #given a set of observations x. Uses the Rcpp function fast.loss.prj
  #to iterate over the observations to get the loss
  K = length(P.list)
  n = dim(x)[1]
  m = dim(x)[2]
  x = as.matrix(x)
  dist = matrix(0,nrow = K, ncol = n)
  for(k in 1:K){
    # dist[k,] = fast.loss.prj(xS = x, PS = P.list[[k]],muS = mu[[k]],nS = n, dS =dim(P.list[[k]])[2],mS =m )
    dist[k,] = fast_loss_prj(n, dim(P.list[[k]])[2], m, P.list[[k]], x, mu[[k]])
  }
  mindist = apply(dist,2,min)
  whichmindist = apply(dist,2,which.min)
  loss = -1*n*temperature *sum(mindist)
  return(list('loss'=loss,'clus'=whichmindist,'dist'=dist))
}
#' @keywords internal
#' @noRd
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#' @keywords internal
#' @noRd
mygibbs.step.c <- function(wi){
  n = nrow(wi)
  k = ncol(wi)
  
  iteration = 1000000
  for (i in 1:iteration){
    labels = rep(0,n)
    for (j in 1:n){
      labels[j] = sample(1:k, 1, prob = as.vector(wi[j,]))
    }
    if (length(unique(labels))==k){
      break
    } 
  }
  return(labels)
}
#' @keywords internal
#' @noRd
mygibbs.step.b <- function(eik, kappa=1){
  n = nrow(eik)
  k = ncol(eik)
  wi = array(0,c(n,k))
  for (i in 1:n){
    tgt = exp(-kappa*as.vector(eik[i,]))
    wi[i,] = tgt/sum(tgt)
  }
  wi[is.na(wi)] = 1
  wi = wi/rowSums(wi)
  return(wi)
}

#' S3 method to predict class label of new data with 'msm' object
#' 
#' Given an instance of \code{msm} class from \code{\link{msm}} function, predict 
#' class label of a new data.
#' 
#' @param object an \code{'msm'} object from \code{\link{msm}} function.
#' @param newdata an \eqn{(m\times p)} matrix of row-stacked observations.
#' @param ... extra parameters (not necessary).
#' 
#' @return a length-\eqn{m} vector of class labels.
#' 
#' @concept utility
#' @export
predict.msm <- function(object, newdata, ...){
  ## PREPARE : EXPLICIT INPUTS
  X = prec_input_matrix(newdata)
  msmlist = list(object)
  
  p1 = base::ncol(X)
  p2 = base::nrow(object$P[[1]])
  
  if (p1!=p2){
    stop("* predict.msm : new data and the 'msm' object have different dimensions.")
  }
  
  output = predict_msm_single(X, msmlist)
  return(as.vector(output))
}

#' @keywords internal
#' @noRd
predict_msm_single <- function(X, msmoutput){
  # modified : list of 'msm' objects 
  # 
  #Given the list of lists of P[[iter]][[k]] as the projection matrices
  #and theta[[iter]][[k]] as the list of thetas
  #assign observations X[n,m] to K clusters. 
  #While we are here, give a confidence interval for the distance of subspaces by calculating the conway distance at each iter
  print.progress = FALSE
  # parameters
  n = nrow(X)
  m = ncol(X)
  iterations = length(msmoutput);
  K = length(msmoutput[[1]]$P)
  
  #store cluster assignments
  cluster.mat = matrix(0,nrow=iterations, ncol=n);
  
  #store distances for each iteration
  dist.mat= matrix(0,nrow=K,ncol=n);
  
  #store each pairwaise distance
  n.pairs = K*(K-1)/2
  subdists = matrix(0,nrow = iterations, ncol = n.pairs);
  # print(paste("Iteration ", 0," at ", Sys.time()))
  for(i in 1:iterations){
    if(i%%10==0){
      if (print.progress){
        print(paste('* msmpred : iteration ', i,'/',iterations,' complete at ', Sys.time(),sep=""))  
      }
    }
    subdist.col = 1 #getting pairwise entry for subdist.col 1 first, incrementing from there
    for(j in 1:(K-1)){
      Pij = msmoutput[[i]]$P[[j]]
      for(h in (j+1):K){
        Pih = msmoutput[[i]]$P[[h]]
        subdists[i,subdist.col] = distance(Pij,Pih,subspace=F)
        # print("distance complete...")
        subdist.col=subdist.col+1;
      }
    }
    #Now get the distance for each observation from the subspace
    # print("42 : start")
    for(k in 1:K){
      P.ik = msmoutput[[i]]$P[[k]]
      theta.ik = (msmoutput[[i]]$theta[[k]])
      dist.mat[k,] = gibbs.loss.prj(x=X,P=P.ik,mu=theta.ik,subspace=F)
    }
    # print("42 : finished")
    #now determine which cluster was closest to the point
    
    cluster.mat[i,]=apply(dist.mat,2,which.min)
  }
  
  # returned object
  # return(list("cluster"=cluster.mat,"sub.dist"=subdists))
  predicted = rep(0,n)
  for (i in 1:n){
    tgt    = as.vector(cluster.mat[,i])
    tbtgt  = table(tgt)
    rnames = rownames(tbtgt)
    predicted[i] = round(as.numeric(rnames[as.numeric(which.max(tbtgt))]))
  }
  # apply(apply(cluster.mat, 2, table),2,which.max)
  
  return(predicted)
}