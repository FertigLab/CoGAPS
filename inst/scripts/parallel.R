devtools::load_all()

data(SimpSim)
nIter <- 1000
system.time({
    results <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3, nEquil=nIter, nSample=nIter)
})
