.onLoad <- function(lib, pkg)
{
  ## Load the rjags wrapper ...
  libdir <- "/home/bst/other/ejfertig/GAPS-JAGS-1.0.2/local/lib"
  if (is.null(getOption("jags.libdir"))) {
      options("jags.libdir" = libdir)
  }
  dyn.load(sprintf('%s/libjags.so', libdir)) 
 
  library.dynam("CoGAPS", pkg, lib, local=FALSE)

  ## ... and the modules
  moddir <- "/home/bst/other/ejfertig/GAPS-JAGS-1.0.2/local/lib/JAGS/modules-1.0.2"
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
