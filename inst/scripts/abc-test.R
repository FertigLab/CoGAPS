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
    #P[3, ] <- P[3, ] / sum(P[3, ])
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

# param: theta
# prior: Gamma(1, 0.1)
# proposal: Gamma(theta^2 / delta^2, theta / delta^2)
prior.shape <- 1
prior.rate <- 0.1
delta <- 5

# param: epsilon
# prior: Exponential(1/3)
# proposal: Exponential(1/epsilon_{i-1})
rate <- 1 / 2

# intialize by sampling theta^{(0)} ~ pi(theta)
theta <- rgamma(iters, prior.shape, prior.rate)
epsilon <- rexp(iters, rate)
ep <- numeric(iters)

for (i in 2:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # simulate theta' ~ K(theta|theta^{(t-1)})
    theta.prime <- rgamma(1, theta[i-1] ^ 2 / delta^2, 
                          theta[i-1] / delta^2)
    # simulate epsilon' ~ K(epsilon|epsilon^{(t-1)})
    eps.prime <- rexp(1, 1 / epsilon[i-1])
    # don't allow epsilons less than 0.25
    eps.prime <- max(eps.prime, 0.25)

    # simulate x ~ p(x | theta')
    growth <- logistic.growth(T, 0, 1, theta.prime)
    growth2 <- logistic.growth(T, 0, 1, 3)
    tmp <- c(growth, growth2)
    #tmp <- tmp / sum(tmp)
    P[3, ] <- tmp
    D.prime <- A %*% P

    # calculate rho(S(x), S(y))
    diff <- D - D.prime
    rho <- norm(diff, "2")
    ep[i] <- rho
    
    # calculate acceptance probability
    alpha <- dgamma(theta.prime, prior.shape, prior.rate, log=TRUE) + 
             dexp(eps.prime, rate, log=TRUE) + 
             dgamma(theta[i-1], theta.prime ^ 2 / delta ^2, 
                    theta.prime / delta  ^ 2, log=TRUE) +
             dexp(epsilon[i-1], 1 / eps.prime, log=TRUE) -
             dgamma(theta[i-1], prior.shape, prior.rate, log=TRUE) - 
             dexp(epsilon[i-1], rate, log=TRUE) - 
             dgamma(theta.prime, theta[i-1] ^ 2 / delta ^ 2, 
                    theta[i-1] / delta  ^ 2, log=TRUE) -
             dexp(eps.prime, 1 / epsilon[i-1], log=TRUE)
      
    alpha <- exp(alpha)
    alpha <- min(1, alpha * (rho < eps.prime))
    
    # accept with probability alpha
    accept <- runif(1) < alpha
    theta[i] <- theta.prime * accept + theta[i-1] * (1-accept)
    epsilon[i] <- eps.prime * accept + epsilon[i-1] * (1-accept)
}

thresh <- 3
plot(density(epsilon))
ts.plot(theta[epsilon < thresh])
# just view the theta's where epsilon is a good range
quantile(theta[epsilon < thresh], c(0.025, 0.5, 0.975))
sum(epsilon < thresh)
