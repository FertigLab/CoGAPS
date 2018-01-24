#' Generate Seeds for Multiple Concurrent Runs
#'
#' @param chains number of seeds to generate (number of chains to run)
#' @param seed positive values are kept, negative values will be overwritten
#'  by a seed generated from the current time
#' @return vector of randomly generated seeds
#' @export
#' @examples
#' generateSeeds(chains=2, seed=-1)
generateSeeds <- function(chains=2, seed=-1)
{
    if (chains < 2 || (as.integer(chains) != chains))
    {
        stop("chains must be >= 2 and an integer")
    }

    if (seed < 1)
    {
        secs <- as.numeric(difftime(Sys.time(), paste(Sys.Date(), "00:00"),
            units="secs"))
        secs <- round(secs)
        seeds <- seq_len(chains) * secs
        return(seeds)
    }
    else
    {
        seeds <- seq_len(chains) * seed
        return(seeds)
    }
}
