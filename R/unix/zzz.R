.onLoad <- function(lib, pkg)
{
  ## Load the rjags wrapper ...
  library.dynam("CoGAPS", pkg, lib, local=FALSE)

  ## ... and the modules
  moddir <- "/usr/local/lib/JAGS/modules-1.0.2"
  if (is.null(getOption("jags.moddir"))) {
      options("jags.moddir" = moddir)
  }
  load.module("basemod")
  load.module("gaps")

  .Call("init_jags_console", PACKAGE="CoGAPS")

  ## Set progress bar type
  if (is.null(getOption("jags.pb"))) {
      options("jags.pb"="text")
  }
}
