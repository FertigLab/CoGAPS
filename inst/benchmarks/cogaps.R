library(CoGAPS)
library(microbenchmark)

runCogapsBenchmark <- function(D, S, k, num)
{
    return(microbenchmark(gapsRun(D, S, nFactor=k, messages=FALSE),
        times=num))
}