library(ggplot)

cogaps.utility <- function() {
    data(SimpSim)
    matplot(t(SimpSim.P))
    x <- prcomp(SimpSim.D)
    matplot(x$rotation[, 1:3])
}

logistic.pattern <- function() {
    data(SimpSim)
    P <- SimpSim.P
    A <- SimpSim.A
    D.old <- SimpSim.D
    S.old <- SimpSim.S

    n <- 10
    T <- seq(-5, 5, length.out = n)

    p3.t <- logistic.growth(T, x.0=0, L=0.3, k=2)
    p3.u <- logistic.growth(T, x.0=1, L=0.3, k=1)

    data <- data.frame(t=c(T, T), y=c(p3.t, p3.u), status=rep(c("treated", "untreated"), each=n))

    ggplot(data, aes(x=t, y=y)) +
        geom_point(aes(colour=status)) +
        geom_line(aes(colour=status))

    P[3, ] <- c(p3.t, p3.u)
}

logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
    output <- L / (1 + exp(- k * (x - x.0)))
    return(output)
}
