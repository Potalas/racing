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

if (idx < 6) {
  ##TSP
  scen_file = "~/aco_tests/scenario-tsp.txt"
  alg_name = "tsp"
} else {
  ##QAP
  scen_file = "~/aco_tests/scenario-qap.txt"
  alg_name = "qap"
}

if (idx %% 6 < 3) {
  param_file = "~/aco_tests/parameters_unconditional.txt"
  use_ego = T
  logname <- "EGO"
} else {
  param_file = "~/aco_tests/parameters.txt"
  use_ego = F
  logname <- "default"
}

if (idx %% 3 == 0) {
  budget = 1000
  nr_seeds = 20
} else if (idx %% 3 == 1) {
  budget = 5000
  nr_seeds = 10
} else {
  budget = 10000
  nr_seeds = 5
}

dir.create(paste0("~/aco_tests/execdir/dir_", idx))

log_dir = "/var/scratch/dlvermet/irace_log_aco/"

s <- readScenario(scen_file)
p <- readParameters(param_file)

for (seed in seq(1,nr_seeds)) {
  scen <- s
  scen$seed <- seed
  scen$maxExperiments <- budget
  scen$execDir <- paste0("~/aco_tests/execdir/dir_", idx)
  scen$logFile <- paste0(log_dir, "log_", logname, "_", alg_name, "_budget_", budget, "_seed_", seed)
  scen <- checkScenario(scen)
  irace(scen, p, use_ego, 1)
}
