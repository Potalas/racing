#######################################
## GENERATE CONFIGURATIONS
#######################################

## When called with an unconditional parameter, it
## must return TRUE
conditionsSatisfied <- function (parameters, partialConfiguration, paramName)
{
  condition <- parameters$conditions[[paramName]]
  # If there is no condition, do not waste time evaluating it.
  if (isTRUE(condition)) return(TRUE)

  v <- eval(condition, as.list(partialConfiguration))
  # Return TRUE if TRUE, FALSE if FALSE or NA
  ## FIXME: If we byte-compile the condition, then we should incorporate the
  ## following into the condition directly. See readForbiddenFile.
  v <- !is.na(v) && v
  return(v)
}

new.empty.configuration <- function(parameters)
{
  newConfigurationsColnames <- c(names(parameters$conditions), ".PARENT.")
  return(setNames(as.list(rep(NA, length(newConfigurationsColnames))),
                  newConfigurationsColnames))
}

get.fixed.value <- function(param, parameters)
{
  value <- parameters$domain[[param]][1]
  type <- parameters$types[[param]]
  if (type == "i") {
    return (as.integer(value))
  } else if (type == "c" || type == "o") {
    return (value)
  } else {
    irace.assert (type == "r")
    return (as.double(value))
  }
}

### Uniform sampling for the initial generation
sampleUniform <- function (parameters, nbConfigurations, digits,
                           forbidden = NULL, repair = NULL)
{
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ))
  empty.configuration <- new.empty.configuration(parameters)

  for (idxConfiguration in seq_len(nbConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      configuration <- empty.configuration
      for (p in seq_along(namesParameters)) {
        currentParameter <- namesParameters[p]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          configuration[[p]] <- NA
          next
        }
        # FIXME: We must be careful because parameters$types does not have the
        # same order as namesParameters, because we sample in the order of the
        # conditions.
        currentType <- parameters$types[[currentParameter]]
        if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled          
        } else if (currentType %in% c("i", "r")) {
          newVal <- sample.unif(currentParameter, parameters, currentType, digits)
        } else {
          irace.assert(currentType %in% c("c","o"))
          possibleValues <- parameters$domain[[currentParameter]]
          newVal <- sample(possibleValues, 1)
        }
        configuration[[p]] <- newVal
      }
      configuration <- as.data.frame(configuration, stringsAsFactors=FALSE)
      if (!is.null(repair)) {
        configuration <- repair(configuration, parameters, digits)
      }

      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
        newConfigurations[idxConfiguration,] <- configuration
        break
      }
      forbidden.retries <- forbidden.retries + 1
    }
    if (forbidden.retries >= 100) {
      irace.error("irace tried 100 times to sample from the model a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
  }
  return (newConfigurations)
}

# bo <<- NULL
# 
# initialize_EGO <- function(parameters, n_job = 10) {
#   source_python("candidate_generation.py")
#   bo <<- py_call(modified_BO, parameters, n_job = as.integer(n_job), n_point = as.integer(10), 
#                 max_infill_eval = as.integer(3500), n_restart = as.integer(5))
# }

## Standard, parallel in python ##
# sampleEGO_orig <- function(parameters, configs, logs, nbNewConfigurations) {
#   # if (is.null(bo)) initialize_EGO(parameters)
#   source_python("candidate_generation.py")
#   
#   bo <- py_call(modified_BO, parameters, n_job = as.integer(10), n_point = as.integer(10),
#                 max_infill_eval = as.integer(3500), n_restart = as.integer(5))
#   if (nrow(logs) > 0) {
#     dt_large <- merge(as.data.table(logs), as.data.table(configs), by.x = 'configuration', by.y = '.ID.')
#     setnames(dt_large, "time", "val")
#     bo$add_data(dt_large)
#   }
#   print("Getting new candidates!")
#   # saveRDS(parameters, "Params_running.rds")
#   # py_save_object(r_to_py(parameters), "Params_running.pkl")
#   # saveRDS(dt_large, "dt_large_running.rds")
#   # py_save_object(r_to_py(dt_large), "dt_large_running.pkl")
#   
#   new_cand <- bo$get_candidates(as.integer(nbNewConfigurations))
#   print("Getting was succesfull, adding NA")
#   #To be consistent with table format, add .parent. column
#   new_cand[[".PARENT."]] = NA
#   print("Adding NA parent succesfull!")
#   return(new_cand)
# }

# V1 - standard, parallel in R ##
sampleEGO_orig <- function(parameters, configs, logs, nbNewConfigurations, repair = F, infill = "MGFI", seed = NULL) {
  # if (is.null(bo)) initialize_EGO(parameters)
  # source_python("candidate_generation.py")
  if (is.null(seed)) seed <- sample.int(1e8, 1)
  
  if (py_module_available("bayes_optim")) {
    bo <- import("bayes_optim")
  }
  else {
    bo <- import_from_path("BayesOpt")
  }
  
  param_unconditional <- parameters
  param_unconditional$conditions <- NULL
  parnames <- param_unconditional$names
  
  if (length(configs) > 0) {
    dt_large <- merge(as.data.table(logs), as.data.table(configs), by.x = 'configuration', by.y = '.ID.')
    values <- as.list(dt_large[, 'time'])$time
  }

  if ("irace_BO" %in% py_list_attributes(bo)) {
    sampler <- bo$irace_BO(param_unconditional, max_infill_eval = as.integer(3500), n_restart = as.integer(5),
                          infill = infill)
    setnames(dt_large, "time", "val")
    if (length(configs) > 0) sampler$add_data(dt_large)
    new_cand <- rbindlist(
      mclapply(
        sample.int(1e8, as.integer(nbNewConfigurations)),
        function(x) {
          sampler$get_candidates(as.integer(1), seed = as.integer(x))
        },
        mc.cores = detectCores()))
  }
  else {
    model <- bo$Surrogate$RandomForest(n_estimators= 500)
    space <- bo$SearchSpace$from_dict(param_unconditional, source = "irace")
    sampler <- bo$ParallelBO(search_space = space, obj_fun = NULL, n_point = as.integer(nbNewConfigurations),
                             acquisition_fun = infill, model = model, DoE_size = as.integer(nbNewConfigurations),
                             random_seed = seed, eval_type = "dataframe")
    if (length(configs) > 0) sampler$tell(dt_large[, ..parnames], values)
    new_cand <- rbindlist(
      mclapply(
        sample.int(1e8, as.integer(nbNewConfigurations)),
        function(x) {
          sampler$ask(as.integer(1), seed = as.integer(x))
        },
        mc.cores = detectCores()))
  }

  new_cand <- new_cand[, ..parnames]
  class(new_cand) <- "data.frame"
  new_cand[[".PARENT."]] = NA
  if (repair) new_cand <- repair_conditional(new_cand, parameters)
  return(new_cand)
  
  # bo <- py_call(modified_BO, param_unconditional, n_job = as.integer(1), n_point = as.integer(nbNewConfigurations),
  #               max_infill_eval = as.integer(3500), n_restart = as.integer(5))

  # if (!is.null(logs) && nrow(logs) > 0) {
  #   dt_large <- merge(as.data.table(logs), as.data.table(configs), by.x = 'configuration', by.y = '.ID.')
  #   setnames(dt_large, "time", "val")
  #   bo$add_data(dt_large)
  # }
  # 
  # print("Getting new candidates - mclapply!")

  # time_used <- system.time({new_cand <- rbindlist(
  #   mclapply(
  #     sample.int(1e8, as.integer(nbNewConfigurations)),
  #     function(x) {
  #       bo$get_candidates(as.integer(1), seed = as.integer(x))
  #       },
  #     mc.cores = detectCores()))
  # })

  # if (debuglevel >= 1)
  # print(paste0("Got new candidates (time used:", time_used, ")"))
# 



  # print("Getting new candidates - lapply!")
  #
  # print(system.time({new_cand <- rbindlist(lapply(seq(nbNewConfigurations), function(x) {
  #   bo$get_candidates(as.integer(1))
  # }))}))

  # class(new_cand) <- "data.frame"
  # # print("Getting was succesfull")
  # new_cand[[".PARENT."]] = NA

  # saveRDS(new_cand, "parallel_cand.rds")
  #
  # new_cand <- bo$get_candidates(as.integer(nbNewConfigurations))
  #
  # saveRDS(new_cand, "old_version_cand.rds")
  # if (repair) new_cand <- repair_conditional(new_cand, parameters)
  # return(new_cand)
}

repair_conditional <- function(new_cand, parameters) {
  for (idx in seq(length(parameters$conditions))) {
    condition = parameters$conditions[[idx]]
    variable = names(parameters$conditions)[[idx]]
    if (!isTRUE(condition)) {
      new_cand[!eval(condition), variable] <- NA
    }
  }
  new_cand
}

## V3 -- sample more than needed (both from best surr and acq) and reduce by dist ##
# sampleEGO <- function(parameters, configs, logs, nbNewConfigurations) {
#   # if (is.null(bo)) initialize_EGO(parameters)
#   source_python("candidate_generation.py")
#   
#   bo <- py_call(modified_BO, parameters, n_job = as.integer(1), n_point = as.integer(nbNewConfigurations),
#                 max_infill_eval = as.integer(3500), n_restart = as.integer(5))
#   
#   if (nrow(logs) > 0) {
#     dt_large <- merge(as.data.table(logs), as.data.table(configs), by.x = 'configuration', by.y = '.ID.')
#     setnames(dt_large, "time", "val")
#     bo$add_data(dt_large)
#   }
#   
#   print("Getting new candidates - mclapply!")
#   
#   time_used <- system.time({new_cand <- rbindlist(
#     mclapply(
#       c(-1*sample.int(1e8, ceiling(0.1 * nbNewConfigurations)), sample.int(1e8, ceiling(1.3 * nbNewConfigurations))), 
#       function(x) {
#         bo$get_candidates(as.integer(1), seed = as.integer(x))
#       }, 
#       mc.cores = detectCores()))
#   })
#   
#   # saveRDS(new_cand, "new_cand_running.rds")
#   
#   temp <- apply(as.matrix(new_cand), 2, as.numeric)
#   
#   pdists <- pdist(temp)
#   rownames(pdists) <- rownames(new_cand)
#   
#   # print(rownames(pdists))
#   # print(nrow(pdists))
#   # print("###################")
#   new_cand2 <- new_cand[seq(nbNewConfigurations), ]
#   
#   tryCatch({
#     while (ncol(pdists) > nbNewConfigurations) {
#       #Not very effiecient, but doesnt matter relative to the other executions
#       idx_rm <- which.min(apply(pdists, 2, function(x) {quantile(x, probs = 0.01)[[1]]}))
#       pdists <- pdists[-idx_rm, -idx_rm]
#     }
#     # print(rownames(pdists))
#     # print(nrow(pdists))
#     new_cand2 <- new_cand[as.numeric(rownames(pdists)), ]
#   }, error = function(e) {
#     print(e)
#     print("in except statement")
#   })
#   new_cand <- new_cand2
#   print(paste0("Got new candidates (time used:", time_used, ")"))
#   
#   class(new_cand) <- "data.frame"
#   new_cand[[".PARENT."]] = NA
#   
#   print(new_cand)
#   print("_________________________")
#   
#   return(new_cand)
# }

## V4 -- 2 stage method ##
# sampleEGO <- function(parameters, configs, logs, nbNewConfigurations) {
#   # if (is.null(bo)) initialize_EGO(parameters)
#   source_python("candidate_generation.py")
#   
#   bo <- py_call(modified_BO, parameters, n_job = as.integer(1), n_point = as.integer(nbNewConfigurations),
#                 max_infill_eval = as.integer(3500), n_restart = as.integer(5))
#   
#   if (nrow(logs) > 0) {
#     dt_large <- merge(as.data.table(logs), as.data.table(configs), by.x = 'configuration', by.y = '.ID.')
#     setnames(dt_large, "time", "val")
#     bo$add_data(dt_large)
#   }
#   
#   print("Getting new candidates - mclapply!")
#   
#   N <- detectCores()
#   
#   nr_tot = ceiling(1.4*nbNewConfigurations/N)*N
#   nr_sur <- ceiling(0.1*nr_tot)
#   nr_acq <- nr_tot - nr_sur
#   
#   time_used <- system.time({
#     new_cand <- rbindlist(
#       mclapply(
#         c(-1*sample.int(1e8, nr_sur), sample.int(1e8, nr_acq)), 
#         function(x) {
#           bo$get_candidates(as.integer(1), seed = as.integer(x))
#         }, 
#         mc.cores = N))
#   })
#   
#   
#   
#   temp <- apply(as.matrix(new_cand), 2, as.numeric)
#   
#   pdists <- pdist(temp)
#   rownames(pdists) <- rownames(new_cand)
#   
#   new_cand2 <- new_cand[seq(nbNewConfigurations), ]
#   
#   tryCatch({
#     while (ncol(pdists) > nbNewConfigurations) {
#       #Not very effiecient, but doesnt matter relative to the other executions
#       idx_rm <- which.min(apply(pdists, 2, function(x) {quantile(x, probs = 0.01)[[1]]}))
#       pdists <- pdists[-idx_rm, -idx_rm]
#     }
#     new_cand2 <- new_cand[as.numeric(rownames(pdists)), ]
#   }, error = function(e) {
#     print(e)
#     print("in except statement")
#   })
#   
#   #Arbitrary limit
#   if (nbNewConfigurations > 10) {
#     bo$add_surrogate_data(new_cand2)
#     nr_tot = ceiling(0.25*nbNewConfigurations/N)*N
#     nr_sur <- ceiling(0.1*nr_tot)
#     nr_acq <- nr_tot - nr_sur
#     
#     time_used <- system.time({
#       new_cand3 <- rbindlist(
#         mclapply(
#           c(-1*sample.int(1e8, nr_sur), sample.int(1e8, nr_acq)), 
#           function(x) {
#             bo$get_candidates(as.integer(1), seed = as.integer(x))
#           }, 
#           mc.cores = N))
#     })
#     new_cand4 <- rbind(new_cand2, new_cand3)
#     
#     temp <- apply(as.matrix(new_cand4), 2, as.numeric)
#     
#     pdists <- pdist(temp)
#     rownames(pdists) <- rownames(new_cand4)
#     
#     new_cand2 <- new_cand4[seq(nbNewConfigurations), ]
#     
#     tryCatch({
#       while (ncol(pdists) > nbNewConfigurations) {
#         #Not very effiecient, but doesnt matter relative to the other executions
#         idx_rm <- which.min(apply(pdists, 2, function(x) {quantile(x, probs = 0.01)[[1]]}))
#         pdists <- pdists[-idx_rm, -idx_rm]
#       }
#       new_cand2 <- new_cand4[as.numeric(rownames(pdists)), ]
#     }, error = function(e) {
#       print(e)
#       print("in except statement")
#     })
#   }
#   
#   new_cand <- new_cand2
#   # print(paste0("Got new candidates (time used:", time_used, ")"))
#   
#   class(new_cand) <- "data.frame"
#   new_cand[[".PARENT."]] = NA
#   
#   # print(new_cand)
#   # print("_________________________")
#   
#   return(new_cand)
# }

# To be called the first time before the second race (with indexIter =
# 2) Nb configurations is the number of configurations at the end
# included the elite ones obtained from the previous iteration
sampleModel <- function (parameters, eliteConfigurations, model,
                         nbNewConfigurations, digits, forbidden = NULL,
                         repair = NULL)
{
  if (nbNewConfigurations <= 0) {
    irace.error ("The number of configurations to generate appears to be negative or zero.")
  }
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbNewConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ))
  empty.configuration <- new.empty.configuration(parameters)
  
  for (idxConfiguration in seq_len(nbNewConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      # Choose the elite which will be the parent.
      indexEliteParent <- sample.int (n = nrow(eliteConfigurations), size = 1,
                                      prob = eliteConfigurations[[".WEIGHT."]])
      eliteParent <- eliteConfigurations[indexEliteParent, ]
      idEliteParent <- eliteParent[[".ID."]]
      configuration <- empty.configuration
      configuration[[".PARENT."]] <- idEliteParent
      
      # Sample a value for every parameter of the new configuration.
      for (p in seq_along(namesParameters)) {
        # FIXME: We must be careful because parameters$types does not
        # have the same order as parameters$conditions. Ideally, we
        # should fix this or make it impossible to confuse them.
        currentParameter <- namesParameters[p]
        currentType <- parameters$types[[currentParameter]]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          # Some conditions are unsatisfied.
          # Should be useless, NA is ?always? assigned when matrix created
          newVal <- NA
          
        } else if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled
        } else if (currentType %in% c("i", "r")) {
          mean <- as.numeric(eliteParent[currentParameter])
          # If there is not value we obtain it from the model
          if (is.na(mean)) mean <- model[[currentParameter]][[as.character(idEliteParent)]][2]
          if (is.na(mean)) {
            # The elite parent does not have any value for this parameter,
            # let's sample uniformly.
            newVal <- sample.unif(currentParameter, parameters, currentType, digits)
                                                                                             
          } else {
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]][1]
            newVal <- sample.norm(mean, stdDev, currentParameter, parameters, currentType, digits)
          }
        } else if (currentType == "o") {
          possibleValues <- paramDomain(currentParameter, parameters)  
          value <- eliteParent[currentParameter]
          
          if (is.na(value)) {
            # The elite parent does not have any value for this
            # parameter, let's sample uniformly
            ## FIXME: We should save the last used parameter in the model and use it here.
            newVal <- sample(possibleValues, 1)
          } else {
            # Find the position within the vector of possible
            # values to determine the equivalent integer.
            mean <- match(value, possibleValues) # Return index of value in array
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]]

            # Sample with truncated normal distribution as an integer.
            # See sample.norm() for an explanation.
            newValAsInt <- floor(rtnorm(1, mean + 0.5, stdDev, lower = 1,
                                        upper = length(possibleValues) + 1L))

            # The probability of this happening is very small, but it can happen.
            if (newValAsInt == length(possibleValues) + 1L)
              newValAsInt <- length(possibleValues)
            
            irace.assert(newValAsInt >= 1L && newValAsInt <= length(possibleValues))
            # Get back to categorical values, find the one corresponding to the
            # newVal
            newVal <- possibleValues[newValAsInt]
          } 
        } else if (currentType == "c") {
          # FIXME: Why is idEliteParent character?
          # FIXME: Why the model is <parameter><Parent>? It makes more sense to be <Parent><parameter>.
          probVector <- model[[currentParameter]][[as.character(idEliteParent)]]
          possibleValues <- paramDomain(currentParameter, parameters)
          newVal <- sample(x = possibleValues, size = 1, prob = probVector)
        } else {
          irace.internal.error("Unexpected condition in sampleModel")
        }
        configuration[[p]] <- newVal
      }
      
      configuration <- as.data.frame(configuration, stringsAsFactors = FALSE)
      if (!is.null(repair)) {
        configuration <- repair(configuration, parameters, digits)
      }
      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
        newConfigurations[idxConfiguration,] <- configuration
        break
      }
      forbidden.retries <- forbidden.retries + 1
    }
    if (forbidden.retries >= 100) {
      irace.error("irace tried 100 times to sample from the model a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
  }
  return (newConfigurations)
}

transform.from.log <- function(x, transf, lowerBound, upperBound)
{
  trLower <- attr(transf, "lower") 
  trUpper <- attr(transf, "upper")
  x <- exp(trLower + (trUpper - trLower) * x)
  return(x)
}

transform.to.log <- function(x, transf, lowerBound, upperBound)
{
  trLower <- attr(transf, "lower") 
  trUpper <- attr(transf, "upper")
  return((log(x) - trLower)/(trUpper - trLower))
}
## How to sample integer values?
#
# The problem: If we have an integer with domain [1,3] and we sample a real value
# and round, then there are more chances of getting 2 than 1 or 3:
# [1, 1,5) -> 1
# [1.5, 2,5) -> 2
# [2.5, 3) -> 3
#
# The solution: Sample in [lowerbound, upperbound + 1], that is, [1, 4], then floor():
# [1, 2) -> 1
# [2, 3) -> 2
# [3, 4) -> 3
#
# Why floor() and not trunc()?
# Because trunc(-1.5) -> -1, while floor(-1.5) -> -2, so for a domain [-3,-1]:
#
# [-3, -2) -> -3
# [-2, -1) -> -2
# [-1, 0)  -> -1
#
# Issue 1: We can sample 4 (upperbound + 1). In that case, we return 3.
#
# Issue 2: When sampling from a truncated normal distribution, the extremes are
# not symmetric.
#
# nsamples <- 100000
# table(floor(rtnorm(nsamples, mean=1, sd=1, lower=1,upper=4)))/nsamples
# table(floor(rtnorm(nsamples, mean=3, sd=1, lower=1,upper=4)))/nsamples
#
# To make them symmetric, we translate by 0.5, so that the mean is at the
# actual center of the interval that will produce the same value after
# truncation, e.g., given an integer value of 1, then mean=1.5, which is at the
# center of [1,2).
#
# nsamples <- 100000
# table(floor(rtnorm(nsamples, mean=1.5, sd=1, lower=1,upper=4)))/nsamples
# table(floor(rtnorm(nsamples, mean=3.5, sd=1, lower=1,upper=4)))/nsamples
#
# The above reasoning also works for log-transformed domains, because 
# floor() happens in the original domain, not in the log-transformed one,
# except for the case of log-transformed negative domains, where we have to
# translate by -0.5.
# 
numeric.value.round <- function(type, value, lowerBound, upperBound, digits)
{  
  irace.assert(is.finite(value))
  if (type == "i") {
    value <- floor(value)
    upperBound <- upperBound - 1L # undo the above for the assert
    # The probability of this happening is very small, but it could happen.
    if (value == upperBound + 1L)
      value <- upperBound
  } else
    value <- round(value, digits)

  irace.assert(value >= lowerBound && value <= upperBound)
  return (value)
}

# Sample value for a numerical parameter.
sample.unif <- function(param, parameters, type, digits = NULL)
{
  lowerBound <- paramLowerBound(param, parameters)
  upperBound <- paramUpperBound(param, parameters)
  transf <- parameters$transform[[param]]
  if (type == "i") {
    # +1 for correct rounding before floor()
    upperBound <- 1L + upperBound
  }
  if (transf == "log") {
    value <- runif(1, min = 0, max = 1)
    value <- transform.from.log(value, transf, lowerBound, upperBound)
  } else {
    value <- runif(1, min = lowerBound, max = upperBound)    
  }
  value <- numeric.value.round(type, value, lowerBound, upperBound, digits)
  return(value)
}

sample.norm <- function(mean, sd, param, parameters, type, digits = NULL)
{
  lowerBound <- paramLowerBound(param, parameters)
  upperBound <- paramUpperBound(param, parameters)
  transf <- parameters$transform[[param]]
  if (type == "i") {
    upperBound <- 1L + upperBound
    # Because negative domains are log-transformed to positive domains.
    mean <- mean + 0.5
  }
  
  if (transf == "log") {
    trMean <- transform.to.log(mean, transf, lowerBound, upperBound)
    value <- rtnorm(1, trMean, sd, lower = 0, upper = 1)
    value <- transform.from.log(value, transf, lowerBound, upperBound)
  } else {
    value <- rtnorm(1, mean, sd, lowerBound, upperBound)
  }

  value <- numeric.value.round(type, value, lowerBound, upperBound, digits)
  return(value)
}
