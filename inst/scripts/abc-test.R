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
# proposal: Normal(theta, delta^2)
prior.sd.a <- 10
prior.mean.a <- 0
delta.a <- 3

# B: 
# prior: Normal(5, 10^2) 
# proposal: Normal(theta, delta^2)
prior.sd.b <- 10
prior.mean.b <- 5
delta.b <- 3

# C: 
# prior: Gamma(0.5, 0.5)
# proposal: Gamma(theta^2 / delta^2, theta / delta^2)
prior.shape.c <- 0.5
prior.rate.c <- 0.5
delta.c <- 3

# D: 
# prior: Gamma(1, 1)
# proposal: Gamma(theta^2 / delta^2, theta / delta^2)
prior.shape.d <- 1
prior.rate.d <- 1
delta.d <- 3

# C: 
# prior: Gamma(0.1, 0.1)
# proposal: Gamma(theta^2 / delta^2, theta / delta^2)
prior.shape.e <- 0.1
prior.rate.e <- 0.1
delta.e <- 3

# make copies of P
P.a <- P.b <- P.c <- P.d <- P.e <- P

# intialize by sampling theta^{(0)} ~ pi(theta)
theta.a <- rnorm(iters, prior.mean.a, prior.sd.a)
theta.b <- rnorm(iters, prior.mean.b, prior.sd.b)
theta.c <- rgamma(iters, prior.shape.c, prior.rate.c)
theta.d <- rgamma(iters, prior.shape.d, prior.rate.d)
theta.e <- rgamma(iters, prior.shape.e, prior.rate.e)

for (i in 2:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # 1. simulate theta' ~ K(theta|theta^{(t-1)})
    theta.prime.a <- rnorm(1, theta.a[i-1], delta.a)
    theta.prime.b <- rnorm(1, theta.b[i-1], delta.b)
    theta.prime.c <- rgamma(1, theta.c[i-1] ^ 2 / delta.c^2, 
                            theta.c[i-1] / delta.c^2)
    theta.prime.d <- rgamma(1, theta.d[i-1] ^ 2 / delta.d^2, 
                            theta.d[i-1] / delta.d^2)
    theta.prime.e <- rgamma(1, theta.e[i-1] ^ 2 / delta.e^2, 
                            theta.e[i-1] / delta.e^2)

    # 2. simulate x ~ p(x | theta')
    growth.a <- logistic.growth(T, 0, 1, theta.prime.a)
    growth.b <- logistic.growth(T, 0, 1, theta.prime.b)
    growth.c <- logistic.growth(T, 0, 1, theta.prime.c)
    growth.d <- logistic.growth(T, 0, 1, theta.prime.d)
    growth.e <- logistic.growth(T, 0, 1, theta.prime.e)

    P.prime.a <- P.a
    P.prime.b <- P.b
    P.prime.c <- P.c
    P.prime.d <- P.d
    P.prime.e <- P.e

    P.prime.a[3, 1:10] <- growth.a
    P.prime.b[3, 1:10] <- growth.b
    P.prime.c[3, 1:10] <- growth.c
    P.prime.d[3, 1:10] <- growth.d
    P.prime.e[3, 1:10] <- growth.e

    # 3. If rho(S(x), S(y)) < epsilon
    D.prime.a <- A %*% P.prime.a
    D.prime.b <- A %*% P.prime.b
    D.prime.c <- A %*% P.prime.c
    D.prime.d <- A %*% P.prime.d
    D.prime.e <- A %*% P.prime.e

    diff.a <- D - D.prime.a
    diff.b <- D - D.prime.b
    diff.c <- D - D.prime.c
    diff.d <- D - D.prime.d
    diff.e <- D - D.prime.e

    rho.a <- norm(diff.a, "2")
    rho.b <- norm(diff.b, "2")
    rho.c <- norm(diff.c, "2")
    rho.d <- norm(diff.d, "2")
    rho.e <- norm(diff.e, "2")
    
    rho.a.accept <- rho.a < epsilon
    rho.b.accept <- rho.b < epsilon
    rho.c.accept <- rho.c < epsilon
    rho.d.accept <- rho.d < epsilon
    rho.e.accept <- rho.e < epsilon

    # a. u ~ U(0, 1)
    u.a <- runif(1, 0, 1) / rho.a.accept
    u.b <- runif(1, 0, 1) / rho.b.accept
    u.c <- runif(1, 0, 1) / rho.c.accept
    u.d <- runif(1, 0, 1) / rho.d.accept
    u.e <- runif(1, 0, 1) / rho.e.accept

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
                     dgamma(theta.c[i-1], theta.prime.c ^ 2 / delta.c ^2, 
                            theta.prime.c / delta.c  ^ 2, log=TRUE) -
                     dgamma(theta.prime.c, theta.c[i-1] ^ 2 / delta.c ^ 2, 
                            theta.c[i-1] / delta.c ^ 2, log=TRUE)
    accept.prob.d <- dgamma(theta.prime.d, prior.shape.d, prior.rate.d, 
                            log=TRUE) - 
                     dgamma(theta.d[i-1], prior.shape.d, prior.rate.d, 
                            log=TRUE) +
                     dgamma(theta.d[i-1], theta.prime.d ^ 2 / delta.d ^2, 
                            theta.prime.d / delta.d  ^ 2, log=TRUE) -
                     dgamma(theta.prime.d, theta.d[i-1] ^ 2 / delta.d ^ 2, 
                            theta.d[i-1] / delta.d ^ 2, log=TRUE)
    accept.prob.e <- dgamma(theta.prime.e, prior.shape.e, prior.rate.e, 
                            log=TRUE) - 
                     dgamma(theta.e[i-1], prior.shape.e, prior.rate.e, 
                            log=TRUE) +
                     dgamma(theta.e[i-1], theta.prime.e ^ 2 / delta.e ^2, 
                            theta.prime.e / delta.e  ^ 2, log=TRUE) -
                     dgamma(theta.prime.e, theta.e[i-1] ^ 2 / delta.e ^ 2, 
                            theta.e[i-1] / delta.e ^ 2, log=TRUE)
    
    accept.prob.a <- exp(accept.prob.a)
    accept.prob.b <- exp(accept.prob.b)
    accept.prob.c <- exp(accept.prob.c)
    accept.prob.d <- exp(accept.prob.d)
    accept.prob.e <- exp(accept.prob.e)
    
    # if u < accept.prob 1 else 0
    # also, if rho >= epsilon then 0
    accept.a <- u.a < accept.prob.a
    accept.b <- u.b < accept.prob.b
    accept.c <- u.c < accept.prob.c
    accept.d <- u.d < accept.prob.d
    accept.e <- u.e < accept.prob.e
    
    # 3b, 3c, 4
    theta.a[i] <- theta.prime.a * accept.a + theta.a[i-1] * (1-accept.a)
    theta.b[i] <- theta.prime.b * accept.b + theta.b[i-1] * (1-accept.b)
    theta.c[i] <- theta.prime.c * accept.c + theta.c[i-1] * (1-accept.c)
    theta.d[i] <- theta.prime.d * accept.d + theta.d[i-1] * (1-accept.d)
    theta.e[i] <- theta.prime.e * accept.e + theta.e[i-1] * (1-accept.e)
}

cat("\n")

library(ggplot2)
library(dplyr)
theme_set(theme_classic())

data <- bind_rows(data_frame(x=1:iters, theta=theta.a, prior="Normal(0, 10^2)"),
                  data_frame(x=1:iters, theta=theta.b, prior="Normal(5, 10^2)"),
                  data_frame(x=1:iters, theta=theta.c, prior="Gamma(0.5, 0.5)"),
                  data_frame(x=1:iters, theta=theta.d, prior="Gamma(1, 1)"),
                  data_frame(x=1:iters, theta=theta.c, prior="Gamma(0.1, 0.1)"))

pdf("~/../Downloads/Prior_comparison.pdf")
ggplot(data, aes(x=x, y=theta)) +
  geom_line(aes(colour=prior, linetype=prior),
            alpha=0.5) +
  geom_smooth(aes(colour=prior),
              se=FALSE) +
  geom_hline(yintercept=4) +
  xlab("MCMC Iteration") +
  ylab(expression(theta)) +
  ggtitle("Comparison of priors with A and P fixed\nOne of two growth parameters estimated")
dev.off()