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
devtools::load_all()

set.seed(20)
iters <- 5000

delta <- 2
prior.mean <- 2
prior.sd <- 2
epsilon <- 0.1

thetaTest <- 2
FP <- data.frame(t(P[3, ]))
treatStatus <- rep(0:0, each=ncol(D))
timeRecorded <- T
results <- gapsTransRun(D, S, nFactor=3, thetaTest, FP, 
                        time.of.sample=timeRecorded, 
                        condition=treatStatus, 
                        nEquil=1000, nSample=1000, 
                        numSnapshots=1, nOutR=100,
                        thin=1,
                        prior.mean=2, prior.sd=2,
                        delta=3, epsilon=0.1)

theta <- results$theta
diffs <- theta[1:(length(theta)-1)] - theta[2:length(theta)]
print(mean(abs(diffs) >= 1e-16))
mean(theta[1001:2000])
