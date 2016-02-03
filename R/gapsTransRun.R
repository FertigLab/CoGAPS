 gapsRunTransformation <- function(D, S, nFactor=3, 
                                   # parameters for transformaiton routine
                                   growth.trans="logistic", 
                                   time.of.sample=1:ncol(D), # default
                                   condition=rep(0, ncol(D)),  # treated/untreated (or treat1 treat2, etc, untreated)
                                   nEquil=nBurn, nSample=nIter)
