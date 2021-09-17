## FUNCTIONS FOR 'ECDF' OBJECTS
#  At the ECDF level,
#   * ecdf_auc : area under the CDF (as in 2003 Monti's consensus clustering)
#
#  At the ELIST Level,
#   * elist_check     : list of 'ecdf' objects
#   * elist_fform     : make a function form on a discretized grid
#   * elist_pdist_*   : pairwise distances 
#                 {ks, wass, lp} and selector
#   * elist_converter : returned objects are 'ecdf'; vectors are re-arranged



# * ecdf_auc --------------------------------------------------------------
#' @keywords internal
#' @noRd
ecdf_auc <- function(xf){
  if (!inherits(xf, "ecdf")){
    stop("* ecdf_auc : your input is not a 'ecdf' object.")
  }
  xf_knots = stats::knots(xf) # sorted knot points (actually, data)
  
  m = length(xf_knots)
  A = 0
  for (i in 2:m){
    A = A + (xf_knots[i]-xf_knots[i-1])*xf(xf_knots[i])
  }
  return(A)
}

# * elist_check -----------------------------------------------------------
#' @keywords internal
#' @noRd
elist_check <- function(elist){
  cond1 = (is.list(elist))
  cond2 = all(unlist(lapply(elist, inherits, "ecdf"))==TRUE)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# * elist_fform -----------------------------------------------------------
#' @keywords internal
#' @noRd
elist_fform <- function(elist){
  nlist = length(elist)
  # compute knot points
  allknots = array(0,c(nlist,2))
  for (i in 1:nlist){
    tgt = stats::knots(elist[[i]])
    allknots[i,] = c(min(tgt), max(tgt))
  }
  mint = min(allknots[,1]) - 0.01
  maxt = max(allknots[,2]) + 0.01
  ssize = min((maxt-mint)/1000, 0.001)
  tseq  = seq(mint, maxt, by=ssize)
  # return the list of y values
  outY = list()
  for (i in 1:nlist){
    tgt       = elist[[i]]
    outY[[i]] = tgt(tseq)
  }
  # return the result
  output = list()
  output$tseq = tseq
  output$fval = outY # list of function values
  return(output)
}

# * elist_pdist series ----------------------------------------------------
#   {ks, wass, lp}
#   wass : http://www-users.math.umn.edu/~bobko001/preprints/2016_BL_Order.statistics_Revised.version.pdf
#' @keywords internal
#' @noRd
elist_pdist_ks <- function(elist){
  trflist  = elist_fform(elist)
  flist = trflist$fval
  nlist = length(flist)
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    fi = flist[[i]]
    for (j in (i+1):nlist){
      fj = flist[[j]]
      theval = max(abs(fi-fj))
      output[i,j] <- output[j,i] <- theval[1]
    }
  }
  return(stats::as.dist(output))
}
#' @keywords internal
#' @noRd
elist_pdist_wass <- function(elist, p){
  nlist = length(elist)
  qseq  = base::seq(from=1e-6, to=1-(1e-6), length.out=8128)
  quants = list() # compute quantile functions first
  for (i in 1:nlist){
    quants[[i]] = as.double(stats::quantile(elist[[i]], qseq))
  }
  
  
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    vali = quants[[i]]
    for (j in (i+1):nlist){
      valj = quants[[j]]
      valij = abs(vali-valj)
      
      if (is.infinite(p)){
        output[i,j] <- output[j,i] <- base::max(valij)
      } else {
        theval <- ((integrate_1d(qseq, valij^p))^(1/p))
        output[i,j] <- output[j,i] <- theval
      }
    }
  } 
  
  return(stats::as.dist(output))
}
#' @keywords internal
#' @noRd
elist_pdist_lp <- function(elist, p){
  nlist = length(elist)
  trflist  = elist_fform(elist)
  flist = trflist$fval
  nlist = length(flist)
  output = array(0,c(nlist,nlist))
  if (is.infinite(p)){
    for (i in 1:(nlist-1)){
      fi = flist[[i]]
      for (j in (i+1):nlist){
        fj = flist[[j]]
        output[i,j] <- output[j,i] <- base::max(base::abs(fi-fj))[1]
      }
    } 
  } else {
    for (i in 1:(nlist-1)){
      fi = flist[[i]]
      for (j in (i+1):nlist){
        fj = flist[[j]]
        theval = ((integrate_1d(trflist$tseq, (abs(fi-fj)^p)))^(1/p))
        output[i,j] <- output[j,i] <- theval
      }
    }
  }
  return(stats::as.dist(output))
}
#' @keywords internal
#' @noRd
elist_pdist <- function(elist, method, p){
  if (all(method=="ks")){
    return(elist_pdist_ks(elist))
  } else if (all(method=="wass")){
    return(elist_pdist_wass(elist, p))
  } else if (all(method=="lp")){
    return(elist_pdist_lp(elist, p))
  }
}


# * elist_converter -------------------------------------------------------
#' @keywords internal
#' @noRd
elist_converter <- function(elist, fname="elist"){
  N      = length(elist)
  output = list()
  for (n in 1:N){
    tgt = elist[[n]]
    if (inherits(tgt, "ecdf")){
      output[[n]] = tgt
    } else {
      tgtvec = as.vector(tgt)
      if ((!any(is.infinite(tgtvec)))&&(!any(is.na(tgtvec)))){
        output[[n]] = stats::ecdf(tgtvec)
      } else {
        stop(paste("* ",fname," : ",n,"-th element from 'elist' is neither an 'ecdf' object nor a vector.",sep=""))
      }
    }
  }
  return(output)
}