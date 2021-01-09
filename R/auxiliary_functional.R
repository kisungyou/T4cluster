## AUXILIARY ROUTINES FOR FUNCTIONAL DATA CLUSTERING
#  faux_pdist_lp : compute pairwise distance by 'Lp' metric



# faux_pdist_lp -----------------------------------------------------------
#' @keywords internal
#' @noRd
faux_pdist_lp <- function(finput, p=2){
  # prepare
  rg    = finput$basis$rangeval
  tsize = round((max(rg)-min(rg))/0.1)
  if (tsize < 1000){
    tsize = 1000
  } else if (tsize > 5000){
    tsize = 5000
  } 
  vecx = seq(from=min(rg), to=max(rg), length.out=tsize)
  vecf = fda::eval.fd(vecx, finput)
  
  output = fpp_pdist_lp(vecx, vecf, p)
  return(output)
}