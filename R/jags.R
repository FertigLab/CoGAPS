jags.model <- function(file, data=sys.frame(sys.parent()), inits,
                       n.chains = 1, n.adapt=1000, nchain)
{

    if (missing(file)) {
        stop("Model file name missing")
    }
    if (!file.exists(file)) {
        stop(paste("Model file \"", file, "\" not found", sep=""))
    }
    if (!missing(nchain)) {
        warning("Argument nchain in jags.model is deprecated. Use n.chains.")
        if (missing(n.chains)) {
            n.chains = nchain
        }
    }
    
    p <- .Call("make_console", PACKAGE="CoGAPS") 
    .Call("check_model", p, file, PACKAGE="CoGAPS")

    varnames <- .Call("get_variable_names", p, PACKAGE="CoGAPS")
    if (is.environment(data)) {
        ##Get a list of numeric objects from the supplied environment
        data <- mget(varnames, envir=data, mode="numeric",
                     ifnotfound=list(NULL))
        ##Strip null entries
        data <- data[!sapply(data, is.null)]
    }
    else if (is.list(data)) {
        v <- names(data)
        if (is.null(v)) {
            stop("data must be a named list")
        }
        if (any(nchar(v)==0)) {
            stop("unnamed variables in data list")
        }
        if (any(duplicated(v))) {
            stop("Duplicated names in data list: ",
                 paste(v[duplicated(v)], collapse=" "))
        }
        relevant.variables <- v %in% varnames
        data <- data[relevant.variables]
        unused.variables <- setdiff(v, varnames)
        for (i in seq(along=unused.variables)) {
            warning("Unused variable \"", unused.variables[i], "\" in data")
        }
    }
    else {
        stop("data must be a list or environment")
    }
    
    .Call("compile", p, data, as.integer(n.chains), TRUE, PACKAGE="CoGAPS")

### Setting initial values

    if (!missing(inits)) {

        checkParameters <- function(inits) {
            if(!is.list(inits))
                return (FALSE)

            inames <- names(inits)
            if (is.null(inames) || any(nchar(inames) == 0))
                return (FALSE)

            if (any(duplicated(inames)))
                return (FALSE)
            
            if (any(inames==".RNG.name")) {
                rngname <- inits[[".RNG.name"]]
                if (!is.character(rngname) || length(rngname) != 1)
                    return (FALSE)
                inits[[".RNG.name"]] <- NULL
            }

            if (!all(sapply(inits, is.numeric)))
                return (FALSE)
            
            return (TRUE)
        }
        
        setParameters <- function(inits, chain) {
            if (!is.null(inits[[".RNG.name"]])) {
                .Call("set_rng_name", p, inits[[".RNG.name"]],
                      as.integer(chain), PACKAGE="CoGAPS")
                inits[[".RNG.name"]] <- NULL
            }
            .Call("set_parameters", p, inits, as.integer(chain),
                  PACKAGE="CoGAPS")
        }
        
        init.values <- vector("list", n.chains)
        
        if (is.function(inits)) {
            if (any(names(formals(inits)) == "chain")) {
                for (i in 1:n.chains) {
                    init.values[[i]] <- inits(chain=i)
                }
            }
            else {
                for (i in 1:n.chains) {
                    init.values[[i]] <- inits()
                }
            }
        }
        else if (is.list(inits)) {

            if (checkParameters(inits)) {
                ## Replicate initial values for all chains
                for (i in 1:n.chains) {
                    init.values[[i]] <- inits
                }
            }
            else {
                if (length(inits) != n.chains) {
                    stop("Length mismatch in inits")
                }
                init.values <- inits
            }
        }
            
        for (i in 1:n.chains) {
            if (!checkParameters(init.values[[i]])) {
                stop("Invalid parameters for chain ", i)
            }
            setParameters(init.values[[i]], i)
        }
    }

    .Call("initialize", p, PACKAGE="CoGAPS")

    model.state <- .Call("get_state", p, PACKAGE="CoGAPS")
    model.data <- .Call("get_data", p, PACKAGE="CoGAPS")
    model.code <- readLines(file, warn=FALSE)
    model <- list("ptr" = function() {p},
                  "data" = function() {model.data},
                  "model" = function() {model.code},
                  "state" = function(internal=FALSE)
                  {
                      if(!internal) {
                          for(i in 1:n.chains) {
                              model.state[[i]][[".RNG.state"]] <- NULL
                              model.state[[i]][[".RNG.name"]] <- NULL
                          }
                      }
                      return(model.state)
                  },
                  "nchain" = function()
                  {
                      .Call("get_nchain", p, PACKAGE="CoGAPS")
                  },
                  "iter" = function()
                  {
                      .Call("get_iter", p, PACKAGE="CoGAPS")
                  },
                  "sync" = function() {
                      
                      model.state <<- .Call("get_state", p, PACKAGE="CoGAPS")
                  },
                  "recompile" = function() {
                      ## Clear the console
                      .Call("clear_console", p, PACKAGE="CoGAPS")
                      p <<- .Call("make_console", PACKAGE="CoGAPS")
                      ## Write the model to a temporary file so we can re-read it
                      mf <- tempfile()
                      writeLines(model.code, mf)
                      .Call("check_model", p, mf, PACKAGE="CoGAPS")
                      unlink(mf)
                      ## Re-compile
                      .Call("compile", p, data, n.chains, FALSE, PACKAGE="CoGAPS")
                      ## Re-initialize
                      if (!is.null(model.state)) {
                          if (length(model.state) != n.chains) {
                              stop("Incorrect number of chains in saved state")
                          }
                          for (i in 1:n.chains) {
                              statei <- model.state[[i]]
                              rng <- statei[[".RNG.name"]]
                              if (!is.null(rng)) {
                                  .Call("set_rng_name", p, rng, i, PACKAGE="CoGAPS")
                                  statei[[".RNG.name"]] <- NULL
                              }
                              .Call("set_parameters", p, statei, i, PACKAGE="CoGAPS")
                          }
                          .Call("initialize", p, PACKAGE="CoGAPS")
                          ## Redo adaptation
                          adapting <- .Call("is_adapting", p, PACKAGE="CoGAPS")
                          if(n.adapt > 0 && adapting) {
                              cat("Adapting\n")
                              .Call("update", p, n.adapt, PACKAGE="CoGAPS")
                              if (!.Call("adapt_off", p, PACKAGE="CoGAPS")) {
                                  warning("Adaptation incomplete");
                              }
                          }
                          model.state <<- .Call("get_state", p, PACKAGE="CoGAPS")
                      }
                      invisible(NULL)
                  })
    class(model) <- "jags"

    if (n.adapt > 0) {
        adapt(model, n.adapt)
    }
    return(model)
}

parse.varname <- function(varname) {

  ## Try to parse string of form "a" or "a[n,p:q,r]" where "a" is a
  ## variable name and n,p,q,r are integers

  v <- try(parse(text=varname, n=1), silent=TRUE)
  if (!is.expression(v) || length(v) != 1)
    return(NULL)

  v <- v[[1]]
  if (is.name(v)) {
    ##Full node array requested
    return(list(name=deparse(v)))
  }
  else if (is.call(v) && identical(deparse(v[[1]]), "[") && length(v) > 2) {
    ##Subset requested
    ndim <- length(v) - 2
    lower <- upper <- numeric(ndim)
    if (any(nchar(sapply(v, deparse)) == 0)) {
      ##We have to catch empty indices here or they will cause trouble
      ##below
      return(NULL)
    }
    for (i in 1:ndim) {
      index <- v[[i+2]]
      if (is.numeric(index)) {
        ##Single index
        lower[i] <- upper[i] <- index
      }
      else if (is.call(index) && length(index) == 3 &&
               identical(deparse(index[[1]]), ":") &&
               is.numeric(index[[2]]) && is.numeric(index[[3]]))
        {
          ##Index range
          lower[i] <- index[[2]]
          upper[i] <- index[[3]]
        }
      else return(NULL)
    }
    if (any(upper < lower))
      return (NULL)
    return(list(name = deparse(v[[2]]), lower=lower, upper=upper))
  }
  return(NULL)
}

parse.varnames <- function(varnames)
{
  names <- character(length(varnames))
  lower <- upper <- vector("list", length(varnames))
  for (i in seq(along=varnames)) {
    y <- parse.varname(varnames[i])
    if (is.null(y)) {
      stop(paste("Invalid variable subset", varnames[i]))
    }
    names[i] <- y$name
    if (!is.null(y$lower)) {
      lower[[i]] <- y$lower
    }
    if (!is.null(y$upper)) {
      upper[[i]] <- y$upper
    }
  }
  return(list(names=names, lower=lower, upper=upper))
}


jags.samples <-
  function(model, variable.names, n.iter, thin=1, type="trace", ...)
{
    if (class(model) != "jags")
      stop("Invalid JAGS model")

    if (!is.character(variable.names) || length(variable.names) == 0)
      stop("variable.names must be a character vector")
     
    if (!is.numeric(n.iter) || length(n.iter) != 1 || n.iter <= 0)
      stop("n.iter must be a positive integer")
    if (!is.character(type))
      stop("type must be a character vector")

    pn <- parse.varnames(variable.names)
    .Call("set_monitors", model$ptr(), pn$names, pn$lower, pn$upper,
          as.integer(thin), type, PACKAGE="CoGAPS")
    update(model, n.iter, ...)
    ans <- .Call("get_monitored_values", model$ptr(), type, PACKAGE="CoGAPS")
    for (i in seq(along=variable.names)) {
      .Call("clear_monitor", model$ptr(), pn$names[i], pn$lower[[i]],
            pn$upper[[i]], type, PACKAGE="CoGAPS")
    }
    return(ans)
}

nchain <- function(model)
{
    if (!inherits(model, "jags"))
      stop("Invalid JAGS model object in nchain")
    
    .Call("get_nchain", model$ptr(), PACKAGE="CoGAPS")
}

load.module <- function(name, path, quiet=FALSE)
{
    if (name %in% list.modules()) {
        ## This is a stop-gap measure as JAGS 2.1.0 does allow you
        ## to load the same module twice. This should be fixed in
        ## later versions.
        return(invisible()) #Module already loaded
    }
    
    if (missing(path)) {
        path = getOption("jags.moddir")
        if (is.null(path)) {
            stop("option jags.moddir is not set")
        }
    }
    if (!is.character(path) || length(path) != 1)
        stop("invalid path")
    if (!is.character(name) || length(name) != 1)
        stop("invalid name")

    file <- file.path(path, paste(name, .Platform$dynlib.ext, sep=""))
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }
    if (!isDLLLoaded(file)) {
        ## We must avoid calling dyn.load twice on the same DLL This
        ## may result in the DLL being unloaded and then reloaded,
        ## which will invalidate pointers to the distributions,
        ## functions and factories in the module.
        dyn.load(file)
    }
    ok <- .Call("load_module", name, PACKAGE="CoGAPS")
    if (!ok) {
        stop("module", name, "not found\n", sep=" ")
    }
    else if (!quiet) {
        cat("module", name, "loaded\n", sep=" ")
    }
    invisible()
}

list.modules <- function()
{
    .Call("get_modules", PACKAGE="CoGAPS");
}

isDLLLoaded <- function(file)
{
    dll.list <- getLoadedDLLs()
    for (i in seq(along=dll.list)) {
        if (dll.list[[i]]["path"][1] == file)
            return(TRUE)
    }
    return(FALSE)
}
