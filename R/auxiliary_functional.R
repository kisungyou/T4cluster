## AUXILIARY ROUTINES FOR FUNCTIONAL DATA CLUSTERING
#  faux_pdist_lp : compute pairwise distance by 'Lp' metric



# faux_pdist_lp -----------------------------------------------------------
#' @keywords internal
#' @noRd
faux_pdist_lp <- function(finput, p=2){
  # prepare
  rg   = finput$basis$rangeval
  vecx = seq(from=min(rg), to=max(rg), length.out=1000)
  vecf = fda::eval.fd(vecx, finput)
  
  output = fpp_pdist_lp(vecx, vecf, p)
  return(output)
}