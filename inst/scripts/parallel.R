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
nIter <- 10000
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

library(ggplot2)
theme_set(theme_classic())

x <- a[, 1]
qplot(x, geom="density")
d = rbind(data.frame(t=1:10000, x=a[, 1], chain="1", m="A"),
          data.frame(t=1:10000, x=a[, 2], chain="2", m="A"),
          data.frame(t=1:10000, x=a[, 3], chain="3", m="A"))
e = rbind(data.frame(t=1:10000, x=p[, 1], chain="1", m="P"),
          data.frame(t=1:10000, x=p[, 2], chain="2", m="P"),
          data.frame(t=1:10000, x=p[, 3], chain="3", m="P"))
data <- rbind(d, e)
ggplot(data, aes(x=t, y=x)) + geom_point(aes(colour=chain)) +
    facet_wrap(~m)

library(dplyr)
library(tidyr)
data2 <- data %>%
    group_by(chain, t) %>%
    summarise(a.plus.p=sum(x),
              a.times.p=prod(x),
              sum.sq.op=sum(x^2),
              sq.sum.sq=sqrt(sum(x^2))) %>%
    gather(operation, value, a.plus.p:sq.sum.sq)

ggplot(data2, aes(x=t, y=value)) + 
    geom_point(aes(colour=chain)) +
    facet_wrap(~operation, scales="free_y")

s <- list(a=a, p=p, c=c)

# calculate gelman-rubin
a.Rhat <- round(gelman.rubin(a), 3)
p.Rhat <- round(gelman.rubin(p), 3)
c.Rhat <- round(gelman.rubin(c), 3)

a.Rhat
p.Rhat
c.Rhat

gelman.rubin(a * p)
gelman.rubin(a + p)

message("A atoms: ", a.Rhat, " P atoms ", p.Rhat,  " chiSq ", c.Rhat)
