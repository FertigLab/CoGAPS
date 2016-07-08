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

# ABC - MCMC
# one parameter, larger delta
# zero out third pattern for now
P[3, 1:10] <- 0
T <- seq(-5, 5, len=10)
set.seed(20)
iters <- 5000
epsilon <- 2

prior.sd <- 10
prior.mean <- 0
delta <- 5

# intialize by sampling theta^{(0)} ~ pi(theta)
theta1 <- rnorm(iters, prior.mean, prior.sd)

for (i in 2:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # 1. simulate theta' ~ K(theta|theta^{(t-1)})
    theta1.prime <- rnorm(1, theta1[i-1], delta)

    # 2. simulate x ~ p(x | theta')
    growth1 <- logistic.growth(T, 0, 1, theta1.prime)
    P.prime <- P
    P.prime[3, 1:10] <- growth1

    # 3. If rho(S(x), S(y)) < epsilon
    D.prime <- A %*% P.prime
    diff <- D - D.prime
    rho <- norm(diff, "2")

    if (rho < epsilon) {
        # a. u ~ U(0, 1)
        u1 <- runif(1, 0, 1)

        # b. if u leq pi(theta')/pi*theta^{(t-1)} times 
        #    K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
        #    theta^{(t)} = theta'
        accept.prob1 <- dnorm(theta1.prime, prior.mean, prior.sd, TRUE) - 
                        dnorm(theta1[i-1], prior.mean, prior.sd, TRUE) +
                        dnorm(theta1[i-1], theta1.prime, delta, TRUE) -
                        dnorm(theta1.prime, theta1[i-1], delta, TRUE)
        accept.prob1 <- exp(accept.prob1)

        if (u1 < accept.prob1) {
            theta1[i] <- theta1.prime
        } else {
            # c. otherwise
            theta1[i] <- theta1[i-1]
        }

    } else {
        # 4. otherwise
        theta1[i] <- theta1[i-1]
    }
}

cat("\n")
summary(theta1)
ts.plot(theta1)
