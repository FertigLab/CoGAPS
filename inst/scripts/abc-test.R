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

# reset P values for third pattern
P[3, ] <- 0

T <- seq(-5, 5, len=10)

# simple ABC for treated theta, all else fixed
set.seed(1)
iters <- 1000
diffs1 <- numeric(iters)
proposals1 <- numeric(iters)
diffs2 <- numeric(iters)
proposals2 <- numeric(iters)
accept1 <- numeric(iters)
accept2 <- numeric(iters)

for (i in 1:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # propose theta, construct data
    theta <- abs(rnorm(1, 0, 10))
    growth <- logistic.growth(T, 0, 1, theta)

    # copy over current P (which is just truth)
    P.prime <- P

    # replace treated group in pattern 3 with proposed data
    P.prime[3, 1:10] <- growth

    # now calculate the data D=A x P and the difference D-D'
    D.prime <- A %*% P.prime
    diff <- D - D.prime

    # now calculate the L2-norm
    d <- norm(diff, "2")

    # save proposal and d for inspections
    diffs1[i] <- d
    proposals1[i] <- theta
    
    # same for second growth parameter
    theta <- abs(rnorm(1, 0, 10))
    growth <- logistic.growth(T, 0, 1, theta)
    P.prime <- P
    P.prime[3, 11:20] <- growth
    D.prime <- A %*% P.prime
    diff <- D - D.prime
    d <- norm(diff, "2")
    diffs2[i] <- d
    proposals2[i] <- theta
}

cat("\n")

library(ggplot2)
data <- rbind(data.frame(diffs=diffs1, proposals=proposals1, which="treat"),
              data.frame(diffs=diffs2, proposals=proposals2, which="untreat"))

truth <- data.frame(true=c(4, 3), which=c("treat", "untreat"))

ggplot(data, aes(proposals, diffs)) +
    geom_point(aes(colour=which)) +
    geom_vline(data=truth, 
               aes(xintercept=true, 
                   colour=which),
               linetype=2)

# now accept/reject abc proposals
set.seed(1)
iters <- 10000
diffs <- numeric(iters)
proposals1 <- numeric(iters)
proposals2 <- numeric(iters)
accept1 <- numeric(iters)
accept2 <- numeric(iters)

for (i in 1:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # copy over current P
    P.prime <- P
  
    # propose theta, construct data
    theta1 <- abs(rnorm(1, 0, 10))
    growth1 <- logistic.growth(T, 0, 1, theta1)
    theta2 <- abs(rnorm(1, 0, 10))
    growth2 <- logistic.growth(T, 0, 1, theta2)
    
    # copy in new growth data
    P.prime[3, ] <- c(growth1, growth2)
    
    # calculate total data with proposals
    D.prime <- A %*% P.prime
    diff <- D - D.prime
    
    # calculate l2 norm
    d <- norm(diff, "2")
    diffs[i] <- d
    proposals1[i] <- theta1
    proposals2[i] <- theta2
    
    if (d < 2) {
      # save P
      P <- P.prime
      
      # save proposals
      accept1[i] <- theta1
      accept2[i] <- theta2
    } # else keep old P
}

mean(diffs < 2)

cat("\n")

# now construct dataset of accepted proposals
data <- rbind(data.frame(theta=accept1[accept1 != 0], which="treat"),
              data.frame(theta=accept2[accept2 != 0], which="untreat"))

ggplot(data, aes(x=theta)) +
  geom_density(aes(fill=which),
               alpha=0.25) +
  geom_vline(data=truth, 
             aes(xintercept=true, 
                 colour=which),
             linetype=2) +
  scale_x_continuous(limits=c(1, 6))

data <- rbind(data.frame(x=proposals1, param="Treated Theta"),
              data.frame(x=proposals2, param="Untreated Theta"),
              data.frame(x=diffs, param="Error"))

ggplot(data, aes(x=x)) +
  geom_density(aes(fill=param),
               alpha=0.25) +
  geom_vline(xintercept=c(2, 3, 4),
             linetype=2)

# ABC - MCMC
# zero out third pattern for now
P[3, ] <- 0
set.seed(1)
iters <- 1000
epsilon <- 2

prior.sd <- 10
prior.mean <- 0
delta <- 1

# intialize by sampling theta^{(0)} ~ pi(theta)
theta1 <- rnorm(iters, prior.mean, prior.sd)
theta2 <- rnorm(iters, prior.mean, prior.sd)

for (i in 2:iters) {
    # message bar
    cat(i, "of", iters, "\r")
  
    # 1. simulate theta' ~ K(theta|theta^{(t-1)})
    theta1.prime <- rnorm(1, theta1[i-1], delta)
    theta2.prime <- rnorm(1, theta2[i-1], delta)

    # 2. simulate x ~ p(x | theta')
    growth1 <- logistic.growth(T, 0, 1, theta1.prime)
    growth2 <- logistic.growth(T, 0, 1, theta2.prime)
    P.prime <- P
    P.prime[3, 1:10] <- growth1
    P.prime[3, 11:20] <- growth2

    # 3. If rho(S(x), S(y)) < epsilon
    D.prime <- A %*% P.prime
    diff <- D - D.prime
    rho <- norm(diff, "2")

    if (rho < epsilon) {
        # a. u ~ U(0, 1)
        u1 <- runif(1, 0, 1)
        u2 <- runif(1, 0, 1)

        # b. if u leq pi(theta')/pi*theta^{(t-1)} times 
        #    K(theta^{(t-1)}|theta')/K(theta'|theta^{(t-1)})
        #    theta^{(t)} = theta'
        accept.prob1 <- dnorm(theta.prime1, prior.mean, prior.sd) / 
                        dnorm(theta1[i-1], prior.mean, prior.sd) *
                        dnorm(theta1[i-1], theta.prime1, delta) /
                        dnorm(theta.prime1, theta1[i-1], delta)

        accept.prob2 <- dnorm(theta.prime2, prior.mean, prior.sd) / 
                        dnorm(theta2[i-1], prior.mean, prior.sd) *
                        dnorm(theta2[i-1], theta.prime2, delta) /
                        dnorm(theta.prime2, theta2[i-1], delta)

        if (u1 < accept.prob1) {
            theta1[i] <- theta.prime1
        } else {
            theta1[i] <- theta1[i-1]
        }

        if (u2 < accept.prob2) {
            theta2[i] <- theta.prime2
        } else {
            theta2[i] <- theta2[i-1]
        }
    } else {
        theta1[i] <- theta1[i-1]
        theta2[i] <- theta2[i-1]
    }
}

cat("\n")
