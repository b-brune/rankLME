library(rankLME)
library(tidyverse)
library(glue)


raw_data_no_intercept <- function(
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
  X <- replicate(p, rnorm(n, 0, sd = sd_X))
  
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





classiclme <- function(
    X, y, Z, g, 
    maxit = 10,
    tol = 1e-6,
    intercepts = list(fixed = TRUE, random = TRUE),
    adjust_re = TRUE,
    weighted = FALSE,
    weight_re = weighted,
    weighting_type = "malahanobis",
    use_outlyingness_weights = TRUE,
    leverage_columns = 1:ncol(X), # -> WAS WILL ICH HIER üBERGEBEN?? 
    mean_function = "hodges_lehmann", sd_function = "Qn_corrected",
    control_mean_sd = list(
      mean_function_arguments_fixed = list(),
      mean_function_arguments_random = list(), 
      sd_function_arguments_fixed = list(),
      sd_function_arguments_random = list()),
    mcd = TRUE
) {
  
  stopifnot(
    "X, and Z need to be matrices" = { is.matrix(X) & is.matrix(Z) },
    "X, y, Z and g need to have the same number of observations" = { 
      nrow(X) == length(y) & 
        nrow(X) == nrow(Z) &
        length(g) == length(y)},
    "g needs to contain at least two different groups" = { length(unique(g)) > 1 }
  )
  
  # Sort the variables by group: 
  ord <- order(g)
  
  y <- y[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  X <- X[ord, , drop = FALSE]
  
  g_old <- g
  g <- g[ord]
  
  # Prepare group matching
  n.groups <- length(unique(g))
  group_sizes <- c(table(g))
  group_limits <- cbind(cumsum(group_sizes) - group_sizes + 1, cumsum(group_sizes))
  
  
  # Assign the control parameters! (Hier wäre vielleicht eine control function
  # cool, die alles was nicht übergeben wird auffüllt!)
  list2env(control_mean_sd, envir = environment())
  
    
    #############################################################################
    #
    # Fit the model with random slopes:
    #
    #############################################################################

    ### !!!!!!!!!!!!!!!!!!
    # Match X and Z (for adjustment of random effects)
    
    # Welche Spalte von X korrespondiert zu welcher Spalte in Z? Sollen wir einfach
    # vorschreiben dass die richtig sortiert sind?
    # For now mal schon!
    Z_matching <- as.logical(apply(Z, 2, function(z) sum(apply(X, 2, function(x) all(z == x)))))
    X_matching <- as.logical(apply(X, 2, function(z) sum(apply(Z, 2, function(x) all(z == x)))))
    # Am besten wäre hier ein lookup table
    # DAS FUNKTIONIERT NICHT WENN WIR MEHRERE RANDOM EFFECTS MATCHEN WOLLEN. INSGESAMT IST DAS JA IRGENDWIE BLÖDSINN,
    # AM BESTEN WÄRE HIER ETWAS ZWEIDIMENSIONALES ODER?
    
    
    # Set the model dimensions
    n <- nrow(X)
    p <- ncol(X) 
    k <- ncol(Z)
    
    
    ### Fit iterative rank estimate
    
    # Fit the first regression
    m <- lm(y ~ X - 1)
    
    
    # Coefficients and predicted values
    beta_hat <- coef(m)  
    y_hat <- X %*% beta_hat
    
    # Calculate the residuals
    marg_residuals <- y - y_hat
    
    
    # Predict the random effects
    groupwise_models <- lapply(1:n.groups, function(i) {
      x = group_limits[i, ]
      r = marg_residuals[x[1]:x[2]]
      ZZ = Z[x[1]:x[2], , drop = FALSE]
      coef = solve(t(ZZ) %*% ZZ) %*% t(ZZ) %*% r
      return(coef)
    })
    
    b_hat <- matrix(unlist(groupwise_models), n.groups, k, byrow = TRUE)
    
    
    # Calculate the residuals
    cond_residuals <- c(unlist(apply(cbind(group_limits, b_hat), 1, function(x) {
        marg_residuals[x[1]:x[2]] - Z[x[1]:x[2], ] %*% matrix(x[3:length(x)])
    })))
    
    
    # Estimates for variance components
    sigma_hat <- sd(cond_residuals)
    theta_hat <- apply(b_hat, 2, sd)
    
    
    beta_init <- beta_hat
    
    l <- 1
    
    repeat {
      
      # Store "old" set of coefficient estimates
      
      beta_old <- beta_hat; theta_old <- theta_hat; sigma_old <- sigma_hat; b_hat_old <- b_hat
      
      # Calculate and build the covariance matrix
      Sigma_sq <- matrix(0, n, n)
      
      if (length(theta_hat) > 1) { 
        Sigma_b <- diag(theta_hat^2) 
      } else { Sigma_b <- theta_hat^2 * diag(1) }
      
      
      for (i in 1:nrow(group_limits)) {
        
        Sigma_sq[group_limits[i, 1]:group_limits[i, 2], 
                 group_limits[i, 1]:group_limits[i, 2]] <- 
          rankLME:::matpow(
            sigma_hat^2 * diag(group_sizes[i]) +
              Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE] %*% 
              Sigma_b %*% 
              t(Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE]), -1/2)
      }
      
      # Rescale y and X
      ###
      ### We downweigh all y's that cause unusually large residuals. This would
      ### also happen in case of leverage points; which is not actually desired,
      ### right??!
      ###
      y_star <- Sigma_sq %*% y
      X_star <- Sigma_sq %*% X
      
      
      # Fit the refined model
      m <- lm(y_star ~ X_star - 1)
      
      beta_hat <- m$coef  
      y_hat_tilde <- X %*% beta_hat #{ if (intercepts$fixed) { beta_hat[-1] } else { beta_hat } }
      
      marg_residuals <- y - y_hat_tilde
      
      groupwise_models <- lapply(1:n.groups, function(i) {
        x = group_limits[i, ]
        r = marg_residuals[x[1]:x[2]]
        ZZ = Z[x[1]:x[2], , drop = FALSE]
        coef = solve(t(ZZ) %*% ZZ) %*% t(ZZ) %*% r
        return(coef)     
      })
      
      
      b_hat <- matrix(unlist(groupwise_models), n.groups, k, byrow = TRUE)
      
      cond_residuals <- c(
        unlist(apply(cbind(group_limits, b_hat), 1, function(x) {
            marg_residuals[x[1]:x[2]] - Z[x[1]:x[2], ] %*% matrix(x[3:length(x)])
        })))
      
      # print(beta_hat)
      
      # Estimates for variance components
      sigma_hat <- sd(cond_residuals)
      theta_hat <- apply(b_hat, 2, sd)
      
      
      if (l >= maxit | 
          sum((beta_old - beta_hat)^2) / sum(beta_old^2) < tol & 
          sum((c(sigma_hat, theta_hat) - c(sigma_old, theta_old))^2) / sum(c(sigma_old, theta_old)^2) < tol
      ) break
      
      l <- l + 1
      
    }
    
  
  # Calculate leverage based on hat matrix
  
  results <- list(
    beta = beta_hat,
    beta_init = beta_init,
    sigma = sigma_hat,
    theta = theta_hat
  )
  
  return(results)
}


comp2 = function(n.groups, group_size) {
  res_rank = replicate(100, 
                       { 
                         d <- raw_data_no_intercept(
                           n = n.groups * group_size, 
                           p = 4, 
                           k = 2, 
                           n.groups = n.groups, 
                           group_size = group_size,
                           ranef_distribution = "normal",
                           error_distribution = "normal",
                           sd_b = 0.5, 
                           sd_eps = 1,
                           sd_X = 2, 
                           beta = c(1, 1, 1, 1)
                         )
                         with(d, ranklme(X, y, Z, g, 
                                         intercepts = list(fixed = FALSE, 
                                                           random = FALSE)))
                       }, simplify = FALSE)
  
  
  res_lme = replicate(100, 
                      { 
                        d <- raw_data_no_intercept(
                          n = n.groups * group_size, 
                          p = 4, 
                          k = 2, 
                          n.groups = n.groups, 
                          group_size = group_size,
                          ranef_distribution = "normal",
                          error_distribution = "normal",
                          sd_b = 0.5, 
                          sd_eps = 1,
                          sd_X = 2, 
                          beta = c(1, 1, 1, 1)
                        )
                        with(d, classiclme(X, y, Z, g))
                      }, simplify = FALSE)
  
  
  process_results <- function(res) {
    data.frame(
      beta_0 = unlist(lapply(res, function(x) x$beta[1])),
      beta_1 = unlist(lapply(res, function(x) x$beta[2])),
      beta_2 = unlist(lapply(res, function(x) x$beta[3])),
      beta_3 = unlist(lapply(res, function(x) x$beta[4])),
      sigma = unlist(lapply(res, "[[", "sigma")),
      theta0 = unlist(lapply(res, function(x) x$theta[1])),
      theta1 = unlist(lapply(res, function(x) x$theta[2]))
    )
  }
  
  res = rbind(
    process_results(res_rank) |> mutate(method = "rankbased"),
    process_results(res_lme) |> mutate(method = "lm")
  ) 
  
  plot = res |>
    pivot_longer(-method) |>
    ggplot(aes(name, value, color=method)) +
    geom_boxplot() +
    geom_hline(yintercept=c(0.5, 1), linetype="dotted") + 
    ggtitle(glue("Comparison of lm- and rankbased algorithm, n_i = {group_size}, g = {n.groups}"))
  
  return(list(data = res, plot = plot))
}


set.seed(711)
res55 = comp2(5, 5)
res55$plot

res205 = comp2(20, 5)
res205$plot

res1010 = comp2(10, 10)
res1010$plot

res2020 = comp2(20, 20)
res2020$plot

save(res55, res205, res1010, res2020, file="comparison_algorithm_with_lm_intercept.Rdata")


##################

