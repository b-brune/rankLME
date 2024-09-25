# Reproduce simulation files:

library(batchtools)
library(rankLME)
library(lme4)
library(robustlmm)
library(dplyr)
library(magrittr)
library(data.table)
library(dplyr)

# Create registry
reg <- makeExperimentRegistry("ranklme_experiments",
  packages = c(
    "lme4",
    "robustlmm",
    "rankLME",
    "robustbase",
    "magrittr",
    "parallel"
  )
)

# reg <- loadRegistry("ranklme_reproduction_registry", writeable = TRUE)
# Specifiy the number of CPUs for computation:

# Specify the number of nodes:
# reg$cluster.functions <- makeClusterFunctionsSocket(ncpus = 120, fs.latency = 65)

##
## Define data generating function
##
create_dataset <- function(data,
                           job, # batchtools-specific arguments
                           n, p, k, n.groups, group_size, ranef_distribution,
                           error_distribution, sd_b, sd_eps,
                           sd_X, beta = rep(1, p), # raw data
                           predictors = "gaussian",
                           outlier_type, outlier_proportion, outlier_level, single_outlier_type,
                           outlier_size, leverage_type, y_outlier_type, 
                           ...) {
  
  # Draw the "raw" dataset
  d <- raw_data(
    n = n, p = p, k = k, n.groups = n.groups, group_size = group_size,
    predictors = predictors,
    ranef_distribution = ranef_distribution,
    error_distribution = error_distribution,
    sd_b = sd_b, sd_eps = sd_eps,
    sd_X = 2, beta = rep(1, p), ...
  )
  
  # Introduce outliers
  d_outl <- switch(outlier_type,
                   none = d,
                   y_outliers = y_outliers(
                     raw_data = d, outlier_proportion = outlier_proportion,
                     outlier_level = outlier_level,
                     single_outlier_type = single_outlier_type,
                     outlier_size = outlier_size, y_outlier_type = y_outlier_type
                   ),
                   leverage = leverage(
                     raw_data = d, outlier_proportion = outlier_proportion,
                     outlier_level = outlier_level,
                     single_outlier_type = single_outlier_type,
                     outlier_size = outlier_size,
                     leverage_type = leverage_type
                   ),
                   breakdown = breakdown(
                     raw_data = d, outlier_proportion = outlier_proportion,
                     outlier_level = outlier_level,
                     single_outlier_type = single_outlier_type
                   )
  )
  
  return(d_outl)
}

addProblem(name = "data_with_outliers", fun = create_dataset)

##
## Define algorithms
##

rank <- function(data, job, instance, sd_function, weighted, which_weights=NULL, ...) {
  d <- instance
  p <- ncol(d$X)
  k <- ncol(d$Z)
  n <- nrow(d$X)
  
  rank = ranklme(
    X = d$X, y = d$y, Z = d$Z, g = d$g,
    sd_function = sd_function, 
    control_mean_sd = list(mean_function_arguments_fixed = list(), 
                           mean_function_arguments_random = list(), 
                           sd_function_arguments_fixed = list(p = (p + d$n.groups)), 
                           sd_function_arguments_random = list()),
    weighted = weighted,
    weighting_type = which_weights,
    maxit = 20
  )
  
  return(rank)
}

addAlgorithm("rankbased", fun = rank)

koller <- function(data, job, instance, ...) {
  d <- instance
  p <- ncol(d$X)
  k <- ncol(d$Z)
  n <- nrow(d$X)
  
  # Make a data matrix
  dat <- with(
    d,
    data.frame(
      magrittr::set_colnames(
        cbind(y, X, Z, g),
        c(
          "y", paste0("X", 1:p - 1),
          paste0("Z", 1:k - 1), "g"
        )
      )
    )
  )
  
  # Build formulas for the functions that need formulas:
  f1 <- paste("y", "~", paste(colnames(dat)[3:(p + 1)], collapse = " + "))
  f2 <- paste("(", "1 +", paste(colnames(dat)[(p + 3):(p + k + 1)], collapse = "+"), "|| g )", collapse = "")
  
  f <- formula(paste(f1, f2, sep = " + "))
  
  koller <- rlmer(f, data = dat)
  
  koller_results <- list(
    beta = fixef(koller),
    variances = c(sigma(koller), unlist(lapply(VarCorr(koller), function(z) attr(z, "stddev")))),
    rmse_b = sqrt(mean((as.matrix(ranef(koller)$g) - t(d$b))^2)),
    rmse_is = sqrt(mean((fitted(koller) - d$y)^2))
  )
  
  return(koller_results)
}

addAlgorithm("koller", fun = koller)


lme <- function(data, job, instance, ...) {
  d <- instance
  p <- ncol(d$X)
  k <- ncol(d$Z)
  n <- nrow(d$X)
  
  # Make a data matrix
  dat <- with(
    d,
    data.frame(
      magrittr::set_colnames(
        cbind(y, X, Z, g),
        c(
          "y", paste0("X", 1:p - 1),
          paste0("Z", 1:k - 1), "g"
        )
      )
    )
  )
  
  # Build formulas for the functions that need formulas:
  f1 <- paste("y", "~", paste(colnames(dat)[3:(p + 1)], collapse = " + "))
  f2 <- paste("(", "1 +", paste(colnames(dat)[(p + 3):(p + k + 1)], collapse = "+"), "|| g )", collapse = "")
  
  f <- formula(paste(f1, f2, sep = " + "))
  
  lme <- lmer(f, data = dat)
  
  lme_results <- list(
    beta = fixef(lme),
    variances = c(sigma(lme), unlist(lapply(VarCorr(lme), function(z) attr(z, "stddev")))),
    rmse_b = sqrt(mean((as.matrix(ranef(lme)$g) - t(d$b))^2)),
    rmse_is = sqrt(mean((fitted(lme) - d$y)^2))
  )
  
  return(lme_results)
}

addAlgorithm("lme", fun = lme)


# Define all settings we are interested in:

##
## Consistency and efficiency: No outliers
##
setting0 = list(
  data_with_outliers = data.table::CJ(
    p = 4, 
    k = 2, 
    sd_eps = 1, 
    sd_b = 0.5, 
    sd_X = 2,
    outlier_type = "none",
    n.groups = c(5, 10, 20, 30, 40, 50), 
    group_size = c(5, 10, 20, 30, 40, 50),
    error_distribution = "normal",
    ranef_distribution = "normal",
    setting = 0
) |>
  mutate(n = n.groups * group_size)
)

algo_settings_efficiency <- list(
  rankbased = data.table::CJ(
    weighted = c(TRUE, FALSE),
    sd_function = c("Qn_corrected", "uniMCD"),
    which_weights = "malahanobis"
  ),
  koller = data.table(),
  lme = data.table()
)

ids_setting0 = addExperiments(prob.designs = setting0, algo.designs = algo_settings_efficiency, repls=200)


##
## Behavior under outliers
##

design <- list(data_with_outliers = cbind(
  p = 4, k = 2, sd_eps = 1, sd_b = 0.5, sd_X = 2,
  bind_rows(
    
    # y-outliers
    data.table::CJ(
      n = 400,
      n.groups = 20, group_size = 20,
      outlier_type = "y_outliers",
      y_outlier_type = c("mean", "variance"),
      outlier_size = c(0, 1, seq(0, 300, 50)[-1]),
      outlier_level = "single",
      single_outlier_type = c("random", "sequential"),
      error_distribution = "normal",
      ranef_distribution = "normal",
      outlier_proportion = c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.3),
      setting = 1
    )[(y_outlier_type == "variance" & outlier_size > 0) | (y_outlier_type == "mean" & outlier_size != 1)],
    
    # leverage
    data.table::CJ(
      n = 400,
      n.groups = 20, group_size = 20,
      outlier_type = "leverage",
      leverage_type = c("mean", "variance"),
      outlier_size = c(0, 1,seq(2, 20, 2)),
      outlier_level = "single",
      single_outlier_type = c("random", "sequential"),
      error_distribution = "normal",
      ranef_distribution = "normal",
      outlier_proportion = c(0.0125, 0.025, 0.05, 0.1, 0.2, 0.3),
      setting = 2
    )[(leverage_type == "variance" & outlier_size > 0) | (leverage_type == "mean" & outlier_size != 1)],
    
    # y-outliers breakdown
    data.table::CJ(
      n = 400,
      outlier_type = "y_outliers",
      n.groups = 20, group_size = 20,
      outlier_size = 1000,
      outlier_level = "single",
      single_outlier_type = c("random", "sequential"),
      error_distribution = "normal",
      ranef_distribution = "normal",
      outlier_proportion = seq(0, 0.5, 0.05),
      setting = 3
    ),
    
    # breakdown leverage
    data.table::CJ(
      n = 400,
      n.groups = 20, group_size = 20,
      outlier_type = "leverage",
      leverage_type = c("mean", "variance"),
      outlier_size = 100,
      outlier_level = "single",
      single_outlier_type = c("random", "sequential"),
      error_distribution = "normal",
      ranef_distribution = "normal",
      outlier_proportion = c(0, 0.0125, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
      setting = 4
    )
  )
)
)


algo_settings_robustness <- list(
  rankbased = data.table::CJ(
    weighted = c(TRUE, FALSE),
    sd_function = "Qn_corrected",
    which_weights = "malahanobis"
  ),
  koller = data.table(),
  lme = data.table()
)

ids_robustness = addExperiments(prob.designs = design, algo.designs = algo_settings_robustness, repls = 100)


# Test a few jobs:
a = unwrap(getJobTable())
a[(setting == 2 & weighted)]

testJob(129596)



####
####
####

# Save the results once the jobs are done:

all_jobs = unwrap(getJobTable(findDone()))

consistency = ijoin(all_jobs[setting == 0], reduceResultsDataTable(all_jobs[setting == 0]))
save(consistency, file="consistency.Rdata")

youtliers = ijoin(all_jobs[setting == 1], reduceResultsDataTable(all_jobs[setting == 1]))
save(youtliers, file="youtliers.Rdata")

leverage = ijoin(all_jobs[setting == 2], reduceResultsDataTable(all_jobs[setting == 2]))
save(leverage, file="leverage.Rdata")

breakdown = ijoin(all_jobs[setting == 3], reduceResultsDataTable(all_jobs[setting == 3]))
save(breakdown, file="breakdown.Rdata")

breakdown_leverage = ijoin(all_jobs[setting == 4], reduceResultsDataTable(all_jobs[setting == 6]))
save(breakdown_leverage, file="breakdown_leverage.Rdata")
