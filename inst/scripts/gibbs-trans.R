logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
    output <- L / (1 + exp(- k * (x - x.0)))
    return(output)
}

n <- 100
T <- seq(-5, 5, length.out = n)

rate.treat <- 2
rate.untreat <- 1

p3.t <- logistic.growth(T, x.0=0, L=1, k=rate.treat)
p3.u <- logistic.growth(T, x.0=1, L=1, k=rate.untreat)

devtools::load_all("../trans")

y <- c(p3.t, p3.u)
treatStatus <- rep(0:1, each=n)
timeRecorded <- rep(T, 2)

# cross fingers...
testTrans <- test_run(y, treatStatus, timeRecorded, iter=2000)
