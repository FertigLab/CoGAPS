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

# A: 
# prior: Normal(0, 10^2) 
# proposal: Normal(theta, 5^2)
prior.sd.a <- 10
prior.mean.a <- 0
delta.a <- 5

# B: 
# prior: Normal(5, 10^2) 
# proposal: Normal(theta, 5^2)
prior.sd.b <- 10
prior.mean.b <- 5
delta.b <- 5

# C: 
# prior: Gamma(0.5, 0.5)
# proposal: Normal(theta, 5^2)
prior.shape.c <- 0.5
prior.rate.c <- 0.5
delta.c <- 5

# make copies of P
P.a <- P.b <- P.c <- P

# intialize by sampling theta^{(0)} ~ pi(theta)
theta.a <- rnorm(iters, prior.mean.a, prior.sd.a)
theta.b <- rnorm(iters, prior.mean.b, prior.sd.b)
theta.c <- rgamma(iters, prior.shape.c, prior.rate.c)

for (i in 2:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # 1. simulate theta' ~ K(theta|theta^{(t-1)})
    theta.prime.a <- rnorm(1, theta.a[i-1], delta.a)
    theta.prime.b <- rnorm(1, theta.b[i-1], delta.b)
    theta.prime.c <- rnorm(1, theta.c[i-1], delta.c)

    # 2. simulate x ~ p(x | theta')
    growth.a <- logistic.growth(T, 0, 1, theta.prime.a)
    growth.b <- logistic.growth(T, 0, 1, theta.prime.b)
    growth.c <- logistic.growth(T, 0, 1, theta.prime.c)

    P.prime.a <- P.a
    P.prime.b <- P.b
    P.prime.c <- P.c

    P.prime.a[3, 1:10] <- growth.a
    P.prime.b[3, 1:10] <- growth.b
    P.prime.c[3, 1:10] <- growth.c

    # 3. If rho(S(x), S(y)) < epsilon
    D.prime.a <- A %*% P.prime.a
    D.prime.b <- A %*% P.prime.b
    D.prime.c <- A %*% P.prime.c

    diff.a <- D - D.prime.a
    diff.b <- D - D.prime.b
    diff.c <- D - D.prime.c

    rho.a <- norm(diff.a, "2")
    rho.b <- norm(diff.b, "2")
    rho.c <- norm(diff.c, "2")
    
    rho.a.accept <- rho.a < epsilon
    rho.b.accept <- rho.b < epsilon
    rho.c.accept <- rho.c < epsilon

    # a. u ~ U(0, 1)
    u.a <- runif(1, 0, 1) / rho.a.accept
    u.b <- runif(1, 0, 1) / rho.b.accept
    u.c <- runif(1, 0, 1) / rho.c.accept

    # b. if u leq pi(theta')/pi*theta^{(t-1)} times 
    #    K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
    #    theta^{(t)} = theta'
    accept.prob.a <- dnorm(theta.prime.a, prior.mean.a, prior.sd.a, TRUE) - 
                     dnorm(theta.a[i-1], prior.mean.a, prior.sd.a, TRUE) +
                     dnorm(theta.a[i-1], theta.prime.a, delta.a, TRUE) -
                     dnorm(theta.prime.a, theta.a[i-1], delta.a, TRUE)
    accept.prob.b <- dnorm(theta.prime.b, prior.mean.b, prior.sd.b, TRUE) - 
                     dnorm(theta.b[i-1], prior.mean.b, prior.sd.b, TRUE) +
                     dnorm(theta.b[i-1], theta.prime.b, delta.b, TRUE) -
                     dnorm(theta.prime.b, theta.b[i-1], delta.b, TRUE)
    accept.prob.c <- dgamma(theta.prime.c, prior.shape.c, prior.rate.c, 
                            log=TRUE) - 
                     dgamma(theta.c[i-1], prior.shape.c, prior.rate.c, 
                            log=TRUE) +
                     dnorm(theta.c[i-1], theta.prime.c, delta.c, TRUE) -
                     dnorm(theta.prime.c, theta.c[i-1], delta.c, TRUE)
    
    accept.prob.a <- exp(accept.prob.a)
    accept.prob.b <- exp(accept.prob.b)
    accept.prob.c <- exp(accept.prob.c)
    
    # if u < accept.prob 1 else 0
    # also, if rho >= epsilon then 0
    accept.a <- u.a < accept.prob.a
    accept.b <- u.b < accept.prob.b
    accept.c <- u.c < accept.prob.c
    
    # 3b, 3c, 4
    theta.a[i] <- theta.prime.a * accept.a + theta.a[i-1] * (1-accept.a)
    theta.b[i] <- theta.prime.b * accept.b + theta.b[i-1] * (1-accept.b)
    theta.c[i] <- theta.prime.c * accept.c + theta.c[i-1] * (1-accept.c)
}

cat("\n")

par(mfrow=c(3, 1))
ts.plot(theta.a)
ts.plot(theta.b)
ts.plot(theta.c)
