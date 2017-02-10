library(CoGAPS)

# function to simulate logistic growth
logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
    output <- L / (1 + exp(- k * (x - x.0)))
    return(output)
}

# generate P, A, D, and S data
logistic.pattern <- function(rate.treat=2, rate.untreat=1, 
                             normalize=TRUE, noise=TRUE) {
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
    
    if (normalize) {
        Pnorm <- sum(P[3,])
		
        P[3, ] <- P[3, ] / Pnorm
		A[,3] <- A[,3] * Pnorm
    }
    
    M <- A %*% P
    error <- matrix(rnorm(prod(dim(M)), sd=0.1), 
                    nrow=nrow(M),
                    ncol=ncol(M))
    if (noise==TRUE) {
      D <- M + error
      D[D < 0] <- 0
    } else {
      D <- M
    }
    
    S <- matrix(0.1, nrow=nrow(D), ncol=ncol(M))

    out <- list(D=D, S=S, A=A, P=P)
    return(out)
}

distance <- function(theta.prime, A, P, D, normalize=TRUE) {
    n <- 10
    T <- seq(-5, 5, length.out = n)
    if (length(theta.prime) == 1) {
      growth <- logistic.growth(T, 0, 1, theta.prime)
      growth2 <- logistic.growth(T, 0, 1, 3)
    } else {
      growth <- logistic.growth(T, 0, 1, theta.prime[1])
      growth2 <- logistic.growth(T, 0, 1, theta.prime[2])
    }
    
    tmp <- c(growth, growth2)
    
    if (normalize==TRUE) {
        tmp <- tmp / sum(tmp)
    }
    
    P[3, ] <- tmp
    D.prime <- A %*% P
    
    # calculate rho(S(x), S(y))
    diff <- D - D.prime
    rho <- norm(diff, "2")
    
    return(rho)
}

# candidates for theta
candidates <- seq(0, 20, by=0.01)

# make comparison in one variable fixed case

# unnormalized
patts <- logistic.pattern(4, 3, normalize = FALSE, noise = FALSE)
D <- patts$D; P <- patts$P; A <- patts$A
rho.unnorm <- sapply(candidates, distance, A=A, P=P, D=D, normalize=FALSE)

# normalized
patts <- logistic.pattern(4, 3, normalize = TRUE, noise = FALSE)
D <- patts$D; P <- patts$P; A <- patts$A
rho.norm <- sapply(candidates, distance, A=A, P=P, D=D, normalize=TRUE)

# make a plot
library(ggplot2)
library(dplyr)

d <- bind_rows(data_frame(theta=candidates, 
                          l2norm=rho.unnorm, 
                          type="unnormalized"),
               data_frame(theta=candidates,
                          l2norm=rho.norm, 
                          type="normalized"))

pdf("~/../Downloads/l2norm.pdf")
ggplot(d, aes(x=theta, y=l2norm)) +
  geom_line() +
  geom_vline(xintercept=4, colour="grey", linetype=2) +
  facet_wrap(~type, nrow=2) +
  theme_classic()

# make comparison in both variables free

# candidates for theta
candidates <- seq(0, 20, by=0.05)
candidates <- expand.grid(candidates, candidates)

# unnormalized
patts <- logistic.pattern(4, 3, normalize = FALSE, noise = FALSE)
D <- patts$D; P <- patts$P; A <- patts$A
rho.unnorm <- apply(candidates, 1, 
                    distance, A=A, P=P, D=D, normalize=FALSE)

# normalized
patts <- logistic.pattern(4, 3, normalize = TRUE, noise = FALSE)
D <- patts$D; P <- patts$P; A <- patts$A
rho.norm <- apply(candidates, 1, 
                  distance, A=A, P=P, D=D, normalize=TRUE)

# make a plot
library(ggplot2)
library(dplyr)

d <- bind_rows(data_frame(theta.1=candidates[, 1], 
                          theta.2=candidates[, 2],
                          l2norm=rho.unnorm, 
                          type="unnormalized"),
               data_frame(theta.1=candidates[, 1],
                          theta.2=candidates[, 2],
                          l2norm=rho.norm, 
                          type="normalized"))

ggplot(filter(d, type=="unnormalized"),
       aes(theta.1, theta.2)) +
  geom_tile(aes(fill=log(l2norm))) +
  ggtitle("unnormalized") +
  theme_classic()

ggplot(filter(d, type=="normalized"),
       aes(theta.1, theta.2)) +
  geom_tile(aes(fill=log(l2norm))) +
  ggtitle("normalized") +
  theme_classic()

dev.off()
