library(CoGAPS)

# function to simulate logistic growth
logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
    output <- L / (1 + exp(- k * (x - x.0)))
    return(output)
}

# generate P, A, D, and S data
logistic.pattern <- function(rate.treat=2, rate.untreat=1) {
    data(SimpSim)
    P <- SimpSim.P
    A <- SimpSim.A
    D.old <- SimpSim.D
    S.old <- SimpSim.S

    n <- 10
    T <- seq(-5, 5, length.out = n)

    p3.t <- logistic.growth(T, x.0=0, L=1, k=rate.treat)
    p3.u <- logistic.growth(T, x.0=0, L=1, k=rate.untreat)

    P[3, ] <- c(p3.t, p3.u)
    M <- A %*% P
    error <- matrix(rnorm(prod(dim(M)), sd=0.1), 
                    nrow=nrow(M),
                    ncol=ncol(M))
    D <- M + error
    D[D < 0] <- 0
    S <- matrix(0.1, nrow=nrow(D), ncol=ncol(M))

    out <- list(D=D, S=S, A=A, P=P)
    return(out)
}

patts <- logistic.pattern(4, 3)
D <- patts$D
S <- patts$S
P <- patts$P
A <- patts$A

T <- seq(-5, 5, len=10)

# simple ABC for treated theta, all else fixed
set.seed(1)
iters <- 1000
diffs <- numeric(iters)
proposals <- numeric(iters)

for (i in 1:iters) {
    cat(i, "of", iters, "\r")
    # propose theta, construct data
    theta <- abs(rnorm(1, 0, 10))
    growth <- logistic.growth(T, 0, 1, theta)

    # copy over current P (which is just truth)
    P.prime <- P

    # replace treated group in pattern 3 with proposed data
    P.prime[3, 1:10] <- growth

    # now calculate the data D=A x P and the difference D-D'
    D.prime <- A %*% P.prime
    diff <- D - D.prime

    # now calculate the L2-norm
    d <- norm(diff, "2")

    # save proposal and d for inspections
    diffs[i] <- d
    proposals[i] <- theta
}

library(ggplot2)
data <- data.frame(diffs, proposals)
