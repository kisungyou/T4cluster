.pkgenv <- new.env(parent = emptyenv())

.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
  
  # Retrieve Current Version
  this.version = packageVersion("T4cluster")
  
  ## Print on Screen
  packageStartupMessage("** ------------------------------------------------------- **")
  packageStartupMessage("**         T4cluster || Tools for Cluster Analysis          ")
  packageStartupMessage("**")
  packageStartupMessage("** Version    : ",this.version,"       (",this.year,")",sep="")
  packageStartupMessage("** Maintainer : Kisung You  (kisungyou@outlook.com)")
  packageStartupMessage("** Website    : https://kisungyou.com/T4cluster/")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
  packageStartupMessage("** ------------------------------------------------------- **")
}

.onUnload <- function(libpath) {
  library.dynam.unload("T4cluster", libpath)
}