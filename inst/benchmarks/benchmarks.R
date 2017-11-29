# get directory of script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# load packages
library(microbenchmark)
library(CoGAPS)

# load benchmarks
source(paste(sep="/", script.basename, "cogaps.R"))

# run benchmarks
data(SimpSim)
nIter <- 3000
print(runCogapsBenchmark(
    D=SimpSim.D,
    S=SimpSim.S,
    nEquil=nIter,
    nSample=nIter,
    nFactor=7,
    seed=12345,
    reps=20
))