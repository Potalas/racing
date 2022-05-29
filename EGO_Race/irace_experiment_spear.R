#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("At least one argument must be supplied (idx)", call. = FALSE)
} else if (length(args) > 1) {
  stop("At most one argument must be supplied")
}

idx  <- as.integer(args[1])

library(irace)
library(reticulate)
use_python('/home/dlvermet/local/python/3.7.6/bin/python3')


if (idx %% 2 == 0) {
  use_ego = T
  logname <- "EGO"  
} else {
  use_ego = F
  logname <- "default"
}


if (idx < 2) {
  budget = 1000
  nr_seeds = 25
} else {
  budget = 10000
  nr_seeds = 10
} 

print(budget)
print(nr_seeds)
dir.create(paste0("~/spear_test/execdir/dir_", idx))

log_dir = "/var/scratch/dlvermet/irace_log_spear/"

scen_file <- "~/spear_test/scenario.txt"
param_file = "~/spear_test/parameters-mixed.txt"

s <- readScenario(scen_file)
p <- readParameters(param_file)

for (seed in seq(1,nr_seeds)) {
  scen <- s
  scen$seed <- seed
  scen$maxExperiments <- budget
  scen$execDir <- paste0("~/spear_test/execdir/dir_", idx)
  scen$logFile <- paste0(log_dir, "log_", logname, "_budget_", budget, "_seed_", seed)
  scen <- checkScenario(scen)
  irace(scen, p, use_ego, 1)
}
