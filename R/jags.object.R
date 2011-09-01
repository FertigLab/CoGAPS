update.jags <- function(object, n.iter = 1, by, progress.bar, ...)
{
    if (!is.numeric(n.iter) || n.iter < 1) {
        stop("Invalid n.iter")
    }

    if ("update" %in% names(object)) {
        ## Old-style jags.model object
        ## The progress bar was created by the object$update() function
        
        if (missing(by))
            by <- floor(n.iter/50)

        object$update(n.iter, by)
    }
    else {
        ## New jags.model object (in version 1.0.3-6)

        adapting <- .Call("is_adapting", object$ptr(), PACKAGE="CoGAPS") 
        on.exit(object$sync())
        
        if (missing(progress.bar)) {
            progress.bar <- getOption("jags.pb")
        }
        if (!is.null(progress.bar)) {
            match.arg(progress.bar, c("text","gui","none"))
            if (progress.bar=="none")
                progress.bar <- NULL
        }
        
        do.pb <- interactive() && !is.null(progress.bar) && n.iter >= 100
        if (do.pb) {
            start.iter <- object$iter()
            end.iter <- start.iter + n.iter
            pb <- switch(progress.bar,
                         "text" = txtProgressBar(start.iter, end.iter,
                         initial = start.iter, style=3, width=50, 
                         char=ifelse(adapting,"+","*")),
                         "gui" = updatePB(start.iter, end.iter, adapting))
        }
        
        ## Set refresh frequency for progress bar
        if (missing(by) || by <= 0) {
            by <- min(ceiling(n.iter/50), 100)
        }
        else {
            by <- ceiling(by)
        }

        ## Do updates
        n <- n.iter
        while (n > 0) {
            .Call("update", object$ptr(), min(n,by), PACKAGE="CoGAPS")
            if (do.pb) {
                switch(progress.bar,
                       "text" = setTxtProgressBar(pb, object$iter()),
                       "gui" =  setPB(pb, object$iter()))
            }
            n <- n - by
        }
        if (do.pb) {
            close(pb)
        }
    }

    invisible(NULL)
}

adapt <- function(object, n.iter, ...)
{
    if(.Call("is_adapting", object$ptr(), PACKAGE="CoGAPS")) {
        update(object, n.iter, ...)
        if (!.Call("adapt_off", object$ptr(), PACKAGE="CoGAPS")) {
            warning("Adaptation incomplete. Recreate the model with a longer adaptive phase.")
        }
    }
    invisible(NULL)
}
                  
coef.jags <- function(object, chain = 1, ...) {
    if (!is.numeric(chain) || chain < 1 || chain > object$nchain()) {
        stop("Invalid chain")
    }
    object$state(internal=FALSE)[[chain]]
}

variable.names.jags <- function(object, ...) {
    .Call("get_variable_names", object$ptr(), PACKAGE="CoGAPS")
}
