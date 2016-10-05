#' generateSeeds
#'
#' @param chains number of MCMC chains to be used 
#' @param seed numeric indicating whether to generate seed from system clock. Default is -1. 
#'
#' @return vector of randomly generated seeds for use with gapsRun, gapsMapRun, or GWCoGAPS
#' @export
#'
#' @examples \dontrun{
#' generateSeeds(chains=2, seed=-1)
#' }
#' 
#' 
generateSeeds <- function(chains=2, seed=-1){
    if (chains < 2 || (as.integer(chains) != chains)) {
        stop("chains must be >= 2 and an integer")
    }

    if (seed < 1) {
        secs <- as.numeric(difftime(Sys.time(),
                                    paste(Sys.Date(), "00:00"),
                                    units="secs"))
        secs <- round(secs)
        seeds <- seq_len(chains) * secs

        return(seeds)
    } else {
        seeds <- seq_len(chains) * seed

        return(seeds)
    }
}
