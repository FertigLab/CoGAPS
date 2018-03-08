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

rhos <- function(theta, D, A, P, T) {
    # allocate curve (assume same number of recordings for all conditions

    n <- ncol(P)
    P_old_curve = logistic.growth(T, k=theta[1])
    P_prime_curve = logistic.growth(T, k=theta[2])

    Pnorm = sum(P_prime_curve)
    Pnorm_old = sum(P_old_curve)

    P_prime = P
    A_prime = A

    P_prime[3, ] = P_prime_curve / Pnorm
    A_prime[, 3] = A_prime[, 3] * Pnorm / Pnorm_old

    D_prime = A_prime %*% P_prime
    D_orig = A %*% P

    diff = D - D_prime
    rho = norm(diff, type="2")

    orig = D - D_orig
    rho_thresh = norm(orig, type="2")

    return(c(rho, rho_thresh))
}

#rhos(results[1, 1:2], D, A, P, T)

# ABC - MCMC
# one parameter, larger delta
devtools::load_all()

set.seed(20)
iters <- 5000

delta <- 2
prior.mean <- 2
prior.sd <- 2
epsilon <- 0.1

# start at prior mean

# test with fixed proposal
results <- FixedMatrix(D, A, P, T,
                       iters, delta, epsilon,
                       prior.mean, prior.sd,
                       TRUE)

head(results)
accepted <- results[, 3] < results[, 4]
quantile(results[accepted, 1], c(0.025, 0.5, 0.975))
mean(accepted)

# test with MH proposal
results <- FixedMatrix(D, A, P, T,
                       iters, delta, epsilon,
                       prior.mean, prior.sd,
                       FALSE)

head(results)
accepted <- results[, 3] < results[, 4]
quantile(results[accepted, 1], c(0.025, 0.5, 0.975))
mean(accepted)

