library(CoGAPS)
library(microbenchmark)

runCogapsBenchmark <- function(D, S, nEquil, nSample, nFactor, seed, reps)
{
    return(microbenchmark(gapsRun(D, S,
        nEquil=nEquil,
        nSample=nSample,
        nFactor=nFactor,
        seed=seed,
        messages=FALSE),
    times=reps))
}