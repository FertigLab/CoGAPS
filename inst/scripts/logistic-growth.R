library(ggplot2)
library(CoGAPS)

cogaps.utility <- function() {
    data(SimpSim)
    matplot(t(SimpSim.P))
    x <- prcomp(SimpSim.D)
    matplot(x$rotation[, 1:3])
}

logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
    output <- L / (1 + exp(- k * (x - x.0)))
    return(output)
}

logistic.pattern <- function(rate.treat=2, rate.untreat=1, plot=FALSE) {
    data(SimpSim)
    P <- SimpSim.P
    A <- SimpSim.A
    D.old <- SimpSim.D
    S.old <- SimpSim.S

    n <- 10
    T <- seq(-5, 5, length.out = n)

    p3.t <- logistic.growth(T, x.0=0, L=0.3, k=rate.treat)
    p3.u <- logistic.growth(T, x.0=1, L=0.3, k=rate.untreat)

    if (plot) {
        data <- data.frame(t=c(T, T), y=c(p3.t, p3.u), status=rep(c("treated", "untreated"), each=n))

        ggplot(data, aes(x=t, y=y)) +
            geom_point(aes(colour=status)) +
            geom_line(aes(colour=status))
    }

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

logistic.cogaps <- function(D, S) {
    nIter <- 5000
    nBurn <- 5000
    results <- gapsRun(D, S, nFactor=3, nEquil=nBurn, nSample=nIter)
    P <- results$Pmean
    arrayIdx <- 1:ncol(P)
    matplot(arrayIdx, t(P), type='l', lwd=10)
}
