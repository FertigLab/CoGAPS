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
T <- seq(-5, 5, len=10)
set.seed(20)
iters <- 5000
epsilon <- 2

# prior: Gamma(1, 0.1)
# proposal: Gamma(theta^2 / delta^2, theta / delta^2)
prior.shape <- 1
prior.rate <- 0.1
delta <- 10

# intialize by sampling theta^{(0)} ~ pi(theta)
theta <- rgamma(iters, prior.shape, prior.rate)

for (i in 2:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # 1. simulate theta' ~ K(theta|theta^{(t-1)})
    theta.prime <- rgamma(1, theta[i-1] ^ 2 / delta^2, 
                            theta[i-1] / delta^2)

    # 2. simulate x ~ p(x | theta')
    growth <- logistic.growth(T, 0, 1, theta.prime)

    P[3, 1:10] <- growth

    # 3. If rho(S(x), S(y)) < epsilon
    D.prime <- A %*% P

    diff <- D - D.prime

    rho <- norm(diff, "2")
    
    rho.accept <- rho < epsilon

    # a. u ~ U(0, 1)
    u <- runif(1, 0, 1) / rho.accept

    # b. if u leq pi(theta')/pi*theta^{(t-1)} times 
    #    K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
    #    theta^{(t)} = theta'
    accept.prob <- dgamma(theta.prime, prior.shape, prior.rate, 
                            log=TRUE) - 
                     dgamma(theta[i-1], prior.shape, prior.rate, 
                            log=TRUE) +
                     dgamma(theta[i-1], theta.prime ^ 2 / delta ^2, 
                            theta.prime / delta  ^ 2, log=TRUE) -
                     dgamma(theta.prime, theta[i-1] ^ 2 / delta ^ 2, 
                            theta[i-1] / delta ^ 2, log=TRUE)
    
    accept.prob <- exp(accept.prob)
    
    # if u < accept.prob 1 else 0
    # also, if rho >= epsilon then 0
    accept <- u < accept.prob
    
    # 3b, 3c, 4
    theta[i] <- theta.prime * accept + theta[i-1] * (1-accept)
}

ts.plot(theta)