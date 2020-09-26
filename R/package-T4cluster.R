#' Tools for Cluster Analysis
#'
#' @docType package
#' @name T4cluster
#' @aliases package-T4cluster
#' @noRd
#' @importFrom fda eval.fd create.bspline.basis smooth.basis
#' @import Rdpack
#' @import maotai
#' @importFrom Rdimtools do.pca
#' @importFrom stats as.dist dist
#' @importFrom utils packageVersion getFromNamespace 
#' @importFrom Rcpp evalCpp
#' @useDynLib T4cluster
NULL

# # Clustering the well-known "Canadian temperature" data (Ramsay & Silverman)
# basis <- create.bspline.basis(c(0, 365), nbasis=21, norder=4)
# fdobj <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],basis,
#                       fdnames=list("Day", "Station", "Deg C"))$fd
# fdobj <- smooth.basis(day.5, as.vector(CanadianWeather$dailyAv[,1,"Temperature.C"]),basis)$fd
# 
# data1 <- CanadianWeather$dailyAv[,1,"Temperature.C"]
# data2 <- eval.fd(0:365, fdobj)[,1]
# 
# par(mfrow=c(1,2))
# plot(data1)
# plot(data2)