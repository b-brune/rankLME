#
# Draw raw data
#
#' @export
#'
raw_data <- function(
  n, p, k, n.groups, group_size, predictors = "gaussian", ranef_distribution = "normal",
  error_distribution = "normal", sd_b = 0.5, sd_eps = 1, sd_X = 2, 
  beta = rep(1, p), ...
) {
  
  # Match n and group sizes
  if (length(group_size) == 1) { group_size <- rep(group_size, n.groups) }
  
  # If it doesn't work out, change n
  if (sum(group_size) != n) {
    n <- sum(group_size)
    message("`n` adjusted to match desired number of groups and group sizes.")
  } 
  
  
  # Define group matching
  g <- rep(1:n.groups, times = group_size)
  
  group_limits <- cbind(cumsum(group_size) - group_size + 1, cumsum(group_size))
  
  # Draw random effects (columnwise by group)
  b <- switch(
    ranef_distribution,
    normal = replicate(n.groups, rnorm(k, 0, sd = sd_b)),
    skewed = sd_b * replicate(n.groups, (rchisq(k, df = 2) - 2) / 2),
    t3 = sd_b * replicate(n.groups, rt(k, df = 3) / sqrt(3)),
    poisson = sd_b * replicate(n.groups, (rpois(k, 2) - 2) / sqrt(2)),
    bimodal = 1 / 2 * replicate(n.groups, 
                        if (rbinom(1, 1, 0.5) == 0) { rnorm(k, -1, sd_b) } else { rnorm(k, 1, sd_b) })
  )
  
  # Draw X-matrix and define Z-matrix accordingly (without outliers)
  if (predictors == "gaussian") {
    X <- cbind(1, replicate(p - 1, rnorm(n, 0, sd = sd_X)))
  } else if (predictors == "discrete") {
    X <- cbind(1, replicate(p - 1, rpois(n, lambda = 3)))
  } else if (predictors == "mixed") {
    X <- cbind(1, replicate(ceiling((p-1)/2), rnorm(n, 0, sd = sd_X)), replicate(floor((p-1)/2), rbinom(n, 1, prob=0.5)))
  }
  
  Z <- X[, 1:k, drop = FALSE]

  Z_help <- as.matrix(Matrix::bdiag(
    lapply(1:nrow(group_limits), function(i) {
      Z[group_limits[i, 1]:group_limits[i, 2], , drop = F]
    }))
  )
  
  eps <- switch(
    error_distribution,
    normal = rnorm(n, 0, sd = sd_eps),
    skewed = sd_eps * (rchisq(n, df = 2) - 2) / 2,
    very_skewed = sd_eps * (rchisq(n, df = 1) - 1) / sqrt(2), 
    t3 = sd_eps * rt(n, df = 3) / sqrt(3),
    poisson = sd_eps * (rpois(n, 2) - 2) / sqrt(2))
  
  y <- X %*% beta + c(Z_help %*% c(b)) + eps
  
  return(list(X = X, y = y, Z = Z, b = b, g = g, n.groups = n.groups))
}



#
# Helper function to determine the outlying observations
#
#' @export

determine_outlying <- function(g, n.groups, outlier_proportion, 
                               outlier_level, single_outlier_type = NULL) {
  
  n <- length(g)
  
  if(outlier_proportion > 0) {
    
    outlying <- switch(
      outlier_level,
      groupwise = {
        tmp <- 1:floor(n.groups * outlier_proportion) # calculate number of outlying groups
        which(g %in% tmp) # select observations that belong to the respective groups
      }, 
      single = {
        if (single_outlier_type == "random") { 
          sample(1:n, floor(n * outlier_proportion))
        } else if (single_outlier_type == "sequential") {
          1:floor(n * outlier_proportion) 
        }
      }
    )
    
  } else {
    outlying <- NULL
  }
  
  return(outlying)
}

#
# Now small helper functions to introduce outliers to the given set of data:
#
#' @export
y_outliers <- function(raw_data, outlier_proportion, outlier_level, single_outlier_type,
                       outlier_size, y_outlier_type = "mean", ...) {
  
  outlying <- determine_outlying(
    g = raw_data$g, 
    outlier_proportion = outlier_proportion,
    outlier_level = outlier_level,
    single_outlier_type = single_outlier_type,
    n.groups = raw_data$n.groups
  )
  
  raw_data$y[outlying] <- switch(
    y_outlier_type,
    mean = raw_data$y[outlying] + outlier_size, # keep the variance, just 
                                                # move the group up or down
    variance = raw_data$y[outlying] * outlier_size
  )
  
  return(c(raw_data, outlying = list(outlying)))  

}

#' @export
leverage <- function(raw_data, outlier_proportion, outlier_level, single_outlier_type,
                     outlier_size, outlying_columns = NULL, leverage_type = "mean", ...) {
  
  outlying <- determine_outlying(
    g = raw_data$g, 
    outlier_proportion = outlier_proportion,
    outlier_level = outlier_level,
    single_outlier_type = single_outlier_type,
    n.groups = raw_data$n.groups
  )
  
  k <- ncol(raw_data$Z)
  
  if (is.null(outlying_columns)) outlying_columns <- 2:ncol(raw_data$X)
  
  raw_data$X[outlying, outlying_columns] <- switch(
    leverage_type,
    mean     = raw_data$X[outlying, outlying_columns] + outlier_size,
    variance = raw_data$X[outlying, outlying_columns] * outlier_size
  )
  
  raw_data$Z <- raw_data$X[, 1:k]
  
  return(c(raw_data, outlying = list(outlying)))  
}

#' @export
breakdown <- function(raw_data, outlier_proportion, outlier_level, 
                      single_outlier_type, ...) {
  
  outlying <- determine_outlying(
    g = raw_data$g, 
    outlier_proportion = outlier_proportion,
    outlier_level = outlier_level,
    single_outlier_type = single_outlier_type,
    n.groups = raw_data$n.groups
  )
  
  raw_data$y[outlying] <- 1e6
  
  return(c(raw_data, outlying = list(outlying)))
}


