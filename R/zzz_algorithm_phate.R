#' test using PHATE
#' 
#' 
#' 
#' @export
testPHATE <- function(X, k=2, par_nbd = 5, par_alpha=2){
  dX = stats::dist(X)
  func.phate <- utils::getFromNamespace("hidden_PHATE", "maotai")
  out.phate  <- func.phate(dX, nbdk=par_nbd, alpha=par_alpha)
  myk <- round(k)
  
  sqrtP = sqrt(out.phate$P)
  extD  = stats::dist(sqrtP) # extrinsic distance on sphere
  
  func.medoids <- utils::getFromNamespace("hidden_kmedoids", "maotai")
  output = func.medoids(extD, myk)
  return(output$cluster)
}

# dat = T4cluster::genSMILEY()
# XX = dat$data
# yy = dat$label
# 
# p4 = as.integer(testPHATE(XX, k=3, par_nbd=5, par_alpha=10))
# p5 = as.integer(testPHATE(XX, k=4, par_nbd=5, par_alpha=10))
# p6 = as.integer(testPHATE(XX, k=5, par_nbd=5, par_alpha=10))
# 
# par(mfrow=c(2,2), pty="s")
# plot(XX, col=yy, main="true label", cex=0.9, pch=19)
# plot(XX, col=p4, main="k=3", cex=0.9, pch=19)
# plot(XX, col=p5, main="k=4", cex=0.9, pch=19)
# plot(XX, col=p6, main="k=5", cex=0.9, pch=19)
