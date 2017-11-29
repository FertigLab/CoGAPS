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
print(runCogapsBenchmark(SimpSim.D, SimpSim.S, 3, 20))