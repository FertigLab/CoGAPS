library(CoGAPS)

# function to simulate logistic growth
logistic.growth <- function(x, x.0=0.5, L=1, k=1) {
    output <- L / (1 + exp(- k * (x - x.0)))
    return(output)
}

# generate P, A, D, and S data
logistic.pattern <- function(rate.treat=2, time.half=2) {
    data(SimpSim)
    P <- SimpSim.P
    A <- SimpSim.A
    D.old <- SimpSim.D
    S.old <- SimpSim.S

    n <- ncol(P)
    T <- seq(-1 * time.half, time.half, length.out = n)

    p3 <- logistic.growth(T, k=rate.treat)

    P[3, ] <- p3
    
    # normalize
    Pnorm <- sum(P[3,])
	
    P[3, ] <- P[3, ] / Pnorm
	A[,3] <- A[,3] * Pnorm
    
    M <- A %*% P
    error <- matrix(rnorm(prod(dim(M)), sd=0.1), 
                    nrow=nrow(M),
                    ncol=ncol(M))
    D <- M + error
    D[D < 0] <- 0
    S <- matrix(0.1, nrow=nrow(D), ncol=ncol(M))

    out <- list(D=D, S=S, A=A, P=P, T=T)
    return(out)
}

patts <- logistic.pattern(4, 5)
D <- patts$D
S <- patts$S
P <- patts$P
A <- patts$A
T <- patts$T

# ABC - MCMC
# one parameter, larger delta
set.seed(20)
iters <- 5000

delta <- 2
prior.mean <- 4
prior.sd <- 2
epsilon <- 0.01

devtools::load_all()
results <- FixedMatrix(D, A, P, T,
                       iters, delta, epsilon,
                       prior.mean, prior.sd,
                       FALSE)

mean(results[1, ])
sd(results[1, ])
quantile(results[1, ], c(0.025, 0.5, 0.975))
