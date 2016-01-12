source("logistic-growth.R")

library(R2jags)

x.0 <- 0
L <- 1
k <- 2.5
N <- 1000
#data.time <- seq(-3, 3, length.out = N)
data.time <- runif(N, -3, 3)
data.logit <- logistic.growth(data.time, x.0=x.0, L=L, k=k)
qplot(data.time, data.logit)
qplot(data.time, -log(L / data.logit - 1))

pattern.model <- function() {
    # likelihood
    for (i in 1:N) {
        transformation ~ dnorm(mean.norm, var.norm)
        transformation <- log(L / data.logit[i] - 1)
    } 

    # calculate mean and var
    mean.norm <- k * x.0
    var.norm <- 1000 * k ^ 2

    # non-informative prior
    k ~ dnorm(0, 1000)
}

pattern.data <- c("L", "x.0", "data.logit", "N")

pattern.params <- c("k")

pattern.inits <- function() {
    inits <- list()
    inits$k <- rnorm(1, mean=0, sd=10)
    return(inits)
}

pattern.fit <- jags(data=pattern.data, inits=pattern.inits,
                    parameters.to.save=pattern.params, 
                    model.file=pattern.model, n.chains=3,
                    n.iter=10000)

pattern.fit
