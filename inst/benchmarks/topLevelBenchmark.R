# save file time at beginning to ensure uniqueness
file_name <- paste("cogaps_benchmark_", format(Sys.time(), "%m_%d_%y_%H_%M_%OS"),
    '.RData', sep='')

# load packages
library(CoGAPS)

# display package version
print(packageVersion('CoGAPS'))

# benchmark dimensions

M_dimensions <- c(10, 25, 50, 100, 250, 500, 750, 1000, 2500, 5000)
N_dimensions <- c(10, 25, 50, 100, 250, 500, 750, 1000, 2500, 5000)
nFactor_dimensions <- c(3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100)
nIter_dimensions <- c(1000, 1250, 1500, 1750, 2000, 2500, 3000, 4000, 5000)

M_default <- floor(median(M_dimensions))
N_default <- floor(median(N_dimensions))
nFactor_default <- floor(median(nFactor_dimensions))
nIter_default <- floor(median(nIter_default))

seed_default <- 12345
reps_default <- 5

# benchmark function

data(GIST_TS_20084)
runBenchmark <- function(m, n, k, iter, seed, reps)
{
    set.seed(123)
    test_mat_D <- matrix(sample(as.matrix(GIST.D), m * n, replace=T),
        nrow = m, ncol = n)

    test_mat_S <- matrix(sample(as.matrix(GIST.S), m * n, replace=T),
        nrow = m, ncol = n)

    times <- c()
    for (r in 1:reps)
    {
        start_time <- Sys.time()
        gapsRun(test_mat_D, test_mat_S, nEquil=iter, nSample=iter,
            nFactor=k, seed=seed+r, messages=FALSE, sampleSnapshots=FALSE,
            numSnapshots=iter / 2, nOutR = iter / 2)
        times[r] <- as.numeric(Sys.time() - start_time)
    }
    params <- c(m, n, k, iter, seed, reps)
    secs <- c(min(times), mean(times), median(times), max(times))
    return(c(params, secs))
}

# store sample times and parameters
samples <- NULL

# run M benchmark
for (n in M_dimensions)
{
    print(paste("benchmarking M =", n))
    bm <- runBenchmark(n, N_default, nFactor_default, nIter_default, seed_default,
        reps_default)
    samples <- rbind(samples, bm)
}

# run N benchmark
for (n in N_dimensions)
{
    print(paste("benchmarking N =", n))
    bm <- runBenchmark(M_default, n, nFactor_default, nIter_default, seed_default,
        reps_default)
    samples <- rbind(samples, bm)
}

# run nFactor benchmark
for (n in nFactor_dimensions)
{
    print(paste("benchmarking K =", n))
    bm <- runBenchmark(M_default, N_default, n, nIter_default, seed_default,
        reps_default)
    samples <- rbind(samples, bm)
}

# run nIter benchmark
for (n in nIter_dimensions)
{
    print(paste("benchmarking nIter =", n))
    bm <- runBenchmark(M_default, N_default, nFactor_default, n, seed_default,
        reps_default)
    samples <- rbind(samples, bm)
}

rownames(samples) <- NULL
colnames(samples) <- c('M', 'N', 'K', 'nIter', 'seed', 'reps', 'min', 'mean',
    'median', 'max')
cogaps_benchmark <- data.frame(samples)
cogaps_version <- packageVersion('CoGAPS')

save(cogaps_benchmark, cogaps_version, file=file_name)
