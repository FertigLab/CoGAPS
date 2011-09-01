.onLoad <- function(lib, pkg)
{
### First task is to get installation directory of JAGS

    ## Try environment variable first
    jags.home <- Sys.getenv("JAGS_HOME")
    if (nchar(jags.home)==0) {
        keyname <- "SOFTWARE\\JAGS\\JAGS-2.1.0"
        if (identical(.Platform$r_arch, "x64")) {
            keyname <- paste(keyname,"-x64", sep="")
        }
        ## Look for multi-user installation in registry
        regkey <- try(readRegistry(keyname, hive = "HLM", maxdepth = 1),
                      silent = TRUE)
        if (inherits(regkey, "try-error")) {
            ## Look for single-user installation in registry
            regkey <- try(readRegistry(keyname, hive = "HCU", maxdepth = 1),
                          silent = TRUE)
        }
        if (inherits(regkey, "try-error") || is.null(regkey[["InstallDir"]])) {
            ## Give up
            stop("Failed to locate JAGS 2.1.0 installation.")
        }
        jags.home <- regkey[["InstallDir"]]
    }

    
### Add jags.home to the windows PATH, if not already present

    bindir <- file.path(jags.home, "bin")
    path <- Sys.getenv("PATH")
    split.path <- strsplit(path, .Platform$path.sep)$PATH
    if (!any(split.path == bindir)) {
        path <- paste(bindir, path, sep=.Platform$path.sep)
        Sys.setenv("PATH"=path)
    }
    
### Set the module directory, if the option jags.moddir is not already set
    
    if (is.null(getOption("jags.moddir"))) {
        options("jags.moddir" = file.path(jags.home, "modules"))
    }
    library.dynam("CoGAPS", pkg, lib, local=FALSE)
    load.module("basemod")
    load.module("gaps")
    
    .Call("init_jags_console", PACKAGE="CoGAPS")

### Set progress bar type
    
    if (is.null(getOption("jags.pb"))) {
        options("jags.pb"="text")
    }
}
