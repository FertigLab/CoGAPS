library(CoGAPS)

data(GIST_TS_20084)

nIter <- 3000
gapsRun(D=GIST.D[1:75,], S=GIST.S[1:75,], nEquil=nIter, nSample=nIter,
    nFactor=7, seed=12345)
