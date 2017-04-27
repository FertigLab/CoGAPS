devtools::load_all()

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

    out <- list(D=D, S=S, A=A, P=P)
    return(out)
}

patts <- logistic.pattern(4, 3)
D <- patts$D
S <- patts$S
P <- patts$P
A <- patts$A

n <- 10
treatStatus <- rep(0:1, each=n)
timeRecorded <- rep(seq(-5, 5, len=n), 2)
FP <- data.frame(t(P[3, ]))

devtools::load_all()
out <- gapsTransRun(D, S, nFactor=3, FP, 
                    time.of.sample=timeRecorded, condition=treatStatus, 
                    nEquil=2, nSample=1, numSnapshots=1, nOutR=1)
out$ASnapshots[[1]]
#out <- gapsMapRun(D, S, FP, nFactor=3)
