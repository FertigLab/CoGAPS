cogaps.utility <- function() {
    data(SimpSim)
    matplot(t(SimpSim.P))
    x <- prcomp(SimpSim.D)
    matplot(x$rotation[, 1:3])
}

logistic.pattern <- function() {
}

logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
}
