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
    p3.u <- logistic.growth(T, x.0=1, L=1, k=rate.untreat)

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

# load code
devtools::load_all("../..")

# MCMC parameters
nIter <- 5000
nBurn <- 5000

# get matrices
patts <- logistic.pattern()
D <- patts$D
S <- patts$S
P.true <- patts$P
A.true <- patts$A

# set up measurement info
n <- 10
treatStatus <- rep(0:1, each=n)
timeRecorded <- rep(seq(-5, 5, len=n), 2)

#     results <- gapsRunTransformation(D, S, nFactor=3,
#                                      growth.trans="logistic",
#                                      time.of.sample=timeRecorded,
#                                      condition=treatStatus,
#                                      nEquil=nBurn, nSample=nIter)

# testings gapsTransRun
ABins=data.frame()
PBins=data.frame()
# whatever initial guess is, that's your pattern 
# we could inherit 
# FP <- matrix(1, nrow=3, ncol=ncol(D))
nFactor <- 3
simulation_id="simulation"
nEquil = 1000
nSample = 1000
nOutR = 1000
output_atomic = FALSE
# fixedMatrix = "P"
fixedBinProbs = FALSE
fixedDomain = "N"
sampleSnapshots = TRUE
numSnapshots = 100
alphaA = 0.01
nMaxA = 100000
max_gibbmass_paraA = 100.0
alphaP = 0.01
nMaxP = 100000
max_gibbmass_paraP = 100.0
# parameters for transformaiton routine
growth.trans="logistic"
time.of.sample=timeRecorded
condition=treatStatus

# pass all settings to C++ within a list
#    if (is.null(P)) {
Config = c(simulation_id, output_atomic, fixedBinProbs, fixedDomain, sampleSnapshots);

ConfigNums = c(nFactor, nEquil, nSample, nOutR, alphaA, nMaxA, max_gibbmass_paraA,
               alphaP, nMaxP, max_gibbmass_paraP, numSnapshots);

#Begin logic error checking code

geneNames = rownames(D);
sampleNames = colnames(D);

# label patterns as Patt N
patternNames = c("0");
for(i in 1:nFactor)
{
  patternNames[i] = paste('Patt', i);
}

# call to C++ Rcpp code
cogapResult = cogapsTrans(D, S, FP, ABins, PBins, 
                          Config, ConfigNums,
                          time.of.sample, condition)
