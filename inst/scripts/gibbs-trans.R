logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
    output <- L / (1 + exp(- k * (x - x.0)))
    return(output)
}

n <- 100
T <- seq(-5, 5, length.out = n)

rate.treat <- 2
rate.untreat <- 1

p3.t <- logistic.growth(T, x.0=0, L=1, k=rate.treat) / 2
p3.u <- logistic.growth(T, x.0=1, L=1, k=rate.untreat) / 3

devtools::load_all("../trans")

y <- c(p3.t, p3.u)
treatStatus <- rep(0:1, each=n)
timeRecorded <- rep(T, 2)

# cross fingers...
testTrans <- test_run(y, treatStatus, timeRecorded, iter=2000)

# now try with varying t_0 and t_f
testing <- function(t0=-5, tf=5, n=100, rate.treat=2, rate.untreat=1) {
    T <- seq(t0, tf, length.out = n)

    p3.t <- logistic.growth(T, x.0=0, L=1, k=rate.treat) / 2
    p3.u <- logistic.growth(T, x.0=1, L=1, k=rate.untreat) / 3

    devtools::load_all("../trans")

    y <- c(p3.t, p3.u)
    treatStatus <- rep(0:1, each=n)
    timeRecorded <- rep(T, 2)

    # cross fingers...
    test_run(y, treatStatus, timeRecorded, iter=2000)
}

# explore output with varying t0 tf
out <- testing()
colMeans(out$beta1chain)
apply(out$beta1chain, 2, sd)
out <- testing(t0=-1, tf=4)
colMeans(out$beta1chain)
apply(out$beta1chain, 2, sd)
