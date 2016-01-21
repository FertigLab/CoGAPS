source("logistic-growth.R")

library(R2jags)

# use different language
# update_transformation
# add stupid checks for function

x.0 <- 2.5
L <- 1
k <- 5
N <- 100
data.time <- seq(-5, 5, length.out = N)
# data.time <- runif(N, -3, 3)
# have to be careful with 0s and 1s
# use machine epsilon
# L can always be forced to 1
# might be interested in x.0
# always know the timepoints at which the timepoint was collected
# how does this correspond to biological timepoint
# have to be flexible with this
# input variable
# for now is fixed
data.logit <- boot::logit(logistic.growth(data.time, x.0=x.0, L=L, k=k))
# qplot(data.time, data.logit)
# qplot(data.time, boot::logit(data.logit))

pattern.model <- function() {
    # likelihood
    for (i in 1:N) {
        data.logit[i] ~ dnorm(y.hat[i], tau)
        y.hat[i] <- a + b * data.time[i]
    } 

    # slope and intercept
    a ~ dnorm(0, 0.1)
    b ~ dnorm(0, 0.1)
    tau <- pow(sigma, -2)
    sigma ~ dunif(0, 1)
}

pattern.data <- c("data.time", "data.logit", "N")

pattern.params <- c("a", "b")

pattern.inits <- function() {
    inits <- list()
    inits$a <- rnorm(1, mean=0, sd=10)
    inits$b <- rnorm(1, mean=0, sd=10)
    return(inits)
}

pattern.fit <- jags(data=pattern.data, inits=pattern.inits,
                    parameters.to.save=pattern.params, 
                    model.file=pattern.model, n.chains=3,
                    n.iter=1000)

pattern.fit
