
#' create sets of correlated paterns for testing
#' @keywords internal
#'
#' @param n number of patterns in each set
#' @param len length of each pattern
#' @param s number of sets
#' @param jit jitter factor
#' 
#' @return list of sets of correlated patterns
makeCorPatternSets <- function(n, len, s, jit = 0.1) {
    motherSet <- lapply(1:n, function(x) runif(len))
    res <- list()
    for (s in seq_len(s)) {
      set <- lapply(1:n, function(x) jitter(motherSet[[x]], factor = jit))
      set <- do.call(cbind, set)
      res[[s]] <- set
    }
    res
}

