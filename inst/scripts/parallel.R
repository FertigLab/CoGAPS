devtools::load_all()

# function to calculate gelman-rubin diagnostic
# param is matrix where each column is an MCMC chain of estimates
gelman.rubin <- function(param) {
    # mcmc information
    n <- nrow(param) # number of iterations
    m <- ncol(param) # number of chains

    # calculate the mean of the means
    theta.bar.bar <- mean(colMeans(param))

    # within chain variance
    W <- mean(apply(param, 2, var))

    # between chain variance
    B <- n / (m - 1) * sum((colMeans(param) - theta.bar.bar) ^ 2)

    # variance of stationary distribution
    theta.var.hat <- (1 - 1 / n) * W + 1 / n * B

    # Potential Scale Reduction Factor (PSRF)
    R.hat <- sqrt(theta.var.hat / W)

    return(R.hat)
}

data(SimpSim)
nIter <- 5000
nBurn <- 30000
chains <- 3 # num of MCMC chains
patterns <- 3
a <- p <- c <- matrix(NA, nrow=nIter, ncol=chains)           # empty matrix to store values from each chain
A.mean <- array(NA, c(nrow(SimpSim.D), patterns, chains))    # empty cube to store the estimates for A
P.mean <- array(NA, c(patterns, ncol(SimpSim.D), chains))    # empty cube to store the estimates for P
for (i in 1:3) {
    results <- gapsRun(SimpSim.D, SimpSim.S, nFactor=3, nEquil=nBurn, nSample=nIter)
    a[, i] <- results$atomsASamp
    p[, i] <- results$atomsPSamp
    c[, i] <- results$chiSqValues[(nBurn+1):(nIter+nBurn)]
    A.mean[, , i] <- results$Amean
    P.mean[, , i] <- results$Pmean
}

# calculate gelman-rubin
a.Rhat <- round(gelman.rubin(a), 3)
p.Rhat <- round(gelman.rubin(p), 3)
c.Rhat <- round(gelman.rubin(c), 3)

message("A atoms: ", a.Rhat, " P atoms ", p.Rhat,  " chiSq ", c.Rhat)
