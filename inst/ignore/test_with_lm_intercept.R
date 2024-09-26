
lm_test <- function(
    X, y, Z, g, 
    intercepts = list(fixed = TRUE, random = TRUE),
    adjust_re = TRUE, maxit=10, tol = 1e-6, control_mean_sd = 
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
  
  #############################################################################
  #
  # Fit the model with random slopes:
  #
  #############################################################################
  
  
  # Remove a possible intercept column from X and Z (in all cases --> 
  # die fügen wir immer von Hand dazu!)
  if (length(unique(X[, 1])) == 1)  X <- X[, -1, drop = FALSE]
  if (length(unique(Z[, 1])) == 1)  Z <- Z[, -1, drop = FALSE]
  
  
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
  
  if (intercepts$fixed & intercepts$random) X_matching <- c(TRUE, X_matching)
  ### !!!!!!!!!!!!!!!!!!!
  
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X) + intercepts$fixed
  k <- ncol(Z) + intercepts$random
  
  
  ### Fit iterative rank estimate
  
  # Fit the first regression
  m <- lm(y ~ X)
  
  
  # Coefficients and predicted values
  beta_hat <- coef(m)  
  y_hat <- { if (intercepts$fixed) { cbind(1, X) } else { X } } %*% beta_hat
  
  
  # Calculate the residuals
  marg_residuals <- y - y_hat
  
  
  # Predict the random effects
  groupwise_models <- lapply(1:n.groups, function(i) {
    x = group_limits[i, ]
    r = marg_residuals[x[1]:x[2]]
    Z = X = Z[x[1]:x[2], , drop = FALSE]
    lm(r ~ Z)
  })
  
  b_hat <- matrix(unlist(lapply(groupwise_models, coef)), n.groups, k, byrow = TRUE)
  
  
  # Calculate the residuals
  cond_residuals <- c(unlist(apply(cbind(group_limits, b_hat), 1, function(x) {
    if (intercepts$random) {
      marg_residuals[x[1]:x[2]] - cbind(1, Z[x[1]:x[2], ]) %*% matrix(x[3:length(x)])
    } else {
      marg_residuals[x[1]:x[2]] - Z[x[1]:x[2], ] %*% matrix(x[3:length(x)])
    }
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
            { if (intercepts$random) {
              cbind(1, Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE]) } else {
                Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE] }} %*% 
            Sigma_b %*% 
            t({ if (intercepts$random) {
              cbind(1, Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE]) } else {
                Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE] }}), -1/2)
      
    }
    
    # Rescale y and X
    ###
    ### We downweigh all y's that cause unusually large residuals. This would
    ### also happen in case of leverage points; which is not actually desired,
    ### right??!
    ###
    y_star <- Sigma_sq %*%  y
    X_star <- Sigma_sq %*%  X
    
    
    # Fit the refined model
    m <- lm(y_star ~ X_star - 1)
    
    beta_hat <- coef(m)
    y_hat_tilde <- X %*% beta_hat #{ if (intercepts$fixed) { beta_hat[-1] } else { beta_hat } }
    
    marg_residuals <- y - y_hat_tilde
    
    if (intercepts$fixed) {
      
      alpha <- mean(marg_residuals)
      
      marg_residuals <- marg_residuals - alpha
      
      y_hat <- y_hat_tilde + alpha
      
      beta_hat <- c(alpha, beta_hat)
      
    }
    
    groupwise_models <- lapply(1:n.groups, function(i) {
      x = group_limits[i, ]
      r = marg_residuals[x[1]:x[2]]
      Z = X = Z[x[1]:x[2], , drop = FALSE]
      lm(r ~ Z)
    }
    )
    
    b_hat <- matrix(unlist(lapply(groupwise_models, coef)), n.groups, k, byrow = TRUE)
    # groupwise_leverage_weights <- unlist(lapply(groupwise_models, "[[", "weights"))
    
    
    cond_residuals <- c(unlist(apply(cbind(group_limits, b_hat), 1, function(x) {
      if (intercepts$random) {
        marg_residuals[x[1]:x[2]] - cbind(1, Z[x[1]:x[2], ]) %*% matrix(x[3:length(x)])
      } else {
        marg_residuals[x[1]:x[2]] - Z[x[1]:x[2], ] %*% matrix(x[3:length(x)])
      }
    })))
    
    sigma_hat <- sd(cond_residuals)
    theta_hat <- apply(b_hat, 2, sd)
    
    if (l >= maxit | 
        sum((beta_old - beta_hat)^2) / sum(beta_old^2) < tol & 
        sum((c(sigma_hat, theta_hat) - c(sigma_old, theta_old))^2) / sum(c(sigma_old, theta_old)^2) < tol
    ) break
    
    l <- l + 1
  } 
  
  
  if (sum(X_matching) > 0 & adjust_re) {
    
    # Calculate the means of the random effects
    b_means <- apply(b_hat, 2, mean)
    
    # Perform adjustment
    b_hat <- {
      if (k == 1) {
        b_hat - b_means
      } else {
        t(apply(b_hat, 1, "-", b_means))
      } 
    }
    
    beta_hat[X_matching] <- beta_hat[X_matching] + b_means
    
  }
  
  rownames(b_hat) <- unique(g)
  
  
  
  
  results <- list(
    beta = beta_hat,
    beta_init = beta_init,
    b = b_hat,
    sigma = sigma_hat,
    theta = theta_hat)
  
  return(results)
}



##################

comp = function(n.groups, group_size) {
  res_rank = replicate(100, 
                       { 
                         d <- raw_data(
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
                           beta = rep(1, 4)
                         )
                         with(d, ranklme(X, y, Z, g, 
                                         intercepts = list(fixed = TRUE, 
                                                           random = TRUE)))
                       }, simplify = FALSE)
  
  
  res_lme = replicate(100, 
                      { 
                        d <- raw_data(
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
                          beta = rep(1, 4)
                        )
                        with(d, lm_test(X, y, Z, g))
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



load("comparison_algorithm_with_lm_intercept.Rdata")


set.seed(711)
res55 = comp(5, 5)
res55$plot + ggtitle("n_i = 5, g = 5")

res205 = comp(20, 5)
res205$plot

res1010 = comp(10, 10)
res1010$plot

res2020 = comp(20, 20)
res2020$plot

pdf("algorithm_with_lm.pdf", width=7, height=5)
(
  (res55$plot + ggtitle("n_i = 5, g = 5")) + 
    (res205$plot + ggtitle("n_i = 5, g = 20")) +
    (res1010$plot + ggtitle("n_i = 10, g = 10")) +
    (res2020$plot + ggtitle("n_i = 20, g = 20"))
) + plot_layout(nrow=2, guides="collect") & 
  plot_annotation(title="Comparison of algorithm with rankbased regression and lm") & ylim(0, 2)
dev.off()


#####

library(rankLME)
library(tidyverse)

theme_set(theme_minimal())


process_results <- function(res) {
  data.frame(
    beta_0 = unlist(lapply(res, function(x) x$beta[1])),
    beta_1 = unlist(lapply(res, function(x) x$beta[2])),
    beta_2 = unlist(lapply(res, function(x) x$beta[3])),
    beta_3 = unlist(lapply(res, function(x) x$beta[4])),
    sigma = unlist(lapply(res, "[[", "sigma")),
    theta0 = unlist(lapply(res, function(x) x$theta[1])),
    theta1 = unlist(lapply(res, function(x) x$theta[2])),
    index=1:100
  )
}

n.groups = 5; group_size = 5


ranklme_test <- function(
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
    mean_function_sd = "hodges_lehmann",
    sd_function_sigma = "Qn_corr2",
    mcd = TRUE,
    control_mean_sd = list(
      mean_function_arguments_fixed = list(),
      mean_function_arguments_random = list(), 
      sd_function_arguments_fixed = list(),
      sd_function_arguments_random = list())
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
  
  
  # Remove a possible intercept column from X and Z (in all cases --> 
  # die fügen wir immer von Hand dazu!)
  if (length(unique(X[, 1])) == 1)  X <- X[, -1, drop = FALSE]
  if (length(unique(Z[, 1])) == 1)  Z <- Z[, -1, drop = FALSE]
  
  
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
  
  if (intercepts$fixed & intercepts$random) X_matching <- c(TRUE, X_matching)
  ### !!!!!!!!!!!!!!!!!!!
  
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X) + intercepts$fixed
  k <- ncol(Z) + intercepts$random
  
  
  ### Fit iterative rank estimate
  
  # Fit the first regression
  m <- ranklm(X = X, y = y, intercept = intercepts$fixed, weighted = weighted, 
              weighting_type = weighting_type,
              leverage_columns = leverage_columns,
              mean_function = mean_function,
              mean_function_arguments = mean_function_arguments_fixed)
  
  
  # Coefficients and predicted values
  beta_hat <- m$coef  
  y_hat <- { if (intercepts$fixed) { cbind(1, X) } else { X } } %*% beta_hat
  
  weights <- m$weights
  
  # Calculate the residuals
  marg_residuals <- y - y_hat
  
  
  # Predict the random effects
  groupwise_models <- lapply(1:n.groups, function(i) {
    x = group_limits[i, ]
    ranklm(y = marg_residuals[x[1]:x[2]], 
           X = Z[x[1]:x[2], , drop = FALSE], 
           intercept = intercepts$random,  
           weighted = weight_re, 
           weighting_type = weighting_type,
           mcd=mcd,
           mean_function = mean_function_sd,
           mean_function_arguments = mean_function_arguments_random)
  })
  
  b_hat <- matrix(unlist(lapply(groupwise_models, "[[", "coef")), n.groups, k, byrow = TRUE)
  
  groupwise_leverage_weights <- lapply(groupwise_models, "[[", "weights")
  
  
  # Calculate the residuals
  cond_residuals <- c(unlist(apply(cbind(group_limits, b_hat), 1, function(x) {
    if (intercepts$random) {
      marg_residuals[x[1]:x[2]] - cbind(1, Z[x[1]:x[2], ]) %*% matrix(x[3:length(x)])
    } else {
      marg_residuals[x[1]:x[2]] - Z[x[1]:x[2], ] %*% matrix(x[3:length(x)])
    }
  })))
  
  
  # Estimates for variance components
  sigma_hat <- do.call(sd_function_sigma, c(list(cond_residuals), sd_function_arguments_fixed))
  theta_hat <- apply(b_hat, 2, function(x) do.call(sd_function, c(list(x), sd_function_arguments_random)))
  
  outlyingness_weights <- { 
    if (use_outlyingness_weights) {
      sapply(2 * sigma_hat / abs(cond_residuals), min, 1)
    } else { rep(1, n) } }
  
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
        matpow(
          sigma_hat^2 * diag(group_sizes[i]) +
            { if (intercepts$random) {
              cbind(1, Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE]) } else {
                Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE] }} %*% 
            Sigma_b %*% 
            t({ if (intercepts$random) {
              cbind(1, Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE]) } else {
                Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE] }}), -1/2)
      
    }
    
    # Rescale y and X
    ###
    ### We downweigh all y's that cause unusually large residuals. This would
    ### also happen in case of leverage points; which is not actually desired,
    ### right??!
    ###
    y_star <- Sigma_sq %*% diag(outlyingness_weights) %*% y
    X_star <- Sigma_sq %*% diag(outlyingness_weights) %*% X
    
    
    # Fit the refined model
    m <- ranklm(
      y = y_star, 
      X = X_star,
      intercept = FALSE, #intercepts$fixed, 
      weighted = weighted,
      weighting_type = weighting_type,
      leverage_columns = leverage_columns,
      mean_function = mean_function,
      mean_function_arguments = mean_function_arguments_fixed)
    
    beta_hat <- m$coef  
    y_hat_tilde <- X %*% beta_hat #{ if (intercepts$fixed) { beta_hat[-1] } else { beta_hat } }
    
    marg_residuals <- y - y_hat_tilde
    
    if (intercepts$fixed) {
      
      alpha <- do.call(mean_function, c(list(marg_residuals), mean_function_arguments_fixed))
      
      marg_residuals <- marg_residuals - alpha
      
      y_hat <- y_hat_tilde + alpha
      
      beta_hat <- c(alpha, beta_hat)
      
    }
    
    groupwise_models <- lapply(1:n.groups, function(i) {
      x = group_limits[i, ]
      ranklm(y = marg_residuals[x[1]:x[2]], 
             X = Z[x[1]:x[2], , drop = FALSE], 
             intercept = intercepts$random,  
             weighted = weight_re, weights = groupwise_leverage_weights[[i]],
             mean_function = mean_function_sd, mcd = mcd,
             mean_function_arguments = mean_function_arguments_random)
    })
    
    
    b_hat <- matrix(unlist(lapply(groupwise_models, "[[", "coef")), n.groups, k, byrow = TRUE)
    # groupwise_leverage_weights <- unlist(lapply(groupwise_models, "[[", "weights"))
    
    
    cond_residuals <- c(unlist(apply(cbind(group_limits, b_hat), 1, function(x) {
      if (intercepts$random) {
        marg_residuals[x[1]:x[2]] - cbind(1, Z[x[1]:x[2], ]) %*% matrix(x[3:length(x)])
      } else {
        marg_residuals[x[1]:x[2]] - Z[x[1]:x[2], ] %*% matrix(x[3:length(x)])
      }
    })))
    
    # print(beta_hat)
    
    sigma_hat <- do.call(sd_function_sigma, c(list(cond_residuals), sd_function_arguments_fixed))
    theta_hat <- apply(b_hat, 2, function(x) do.call(sd_function, c(list(x), sd_function_arguments_random)))
    
    outlyingness_weights <- { 
      if (use_outlyingness_weights) {
        sapply(2 * sigma_hat / abs(cond_residuals), min, 1)
      } else { rep(1, n) } }
    
    
    if (l >= maxit | 
        sum((beta_old - beta_hat)^2) / sum(beta_old^2) < tol & 
        sum((c(sigma_hat, theta_hat) - c(sigma_old, theta_old))^2) / sum(c(sigma_old, theta_old)^2) < tol
    ) break
    
    l <- l + 1
    
  }
  
  if (sum(X_matching) > 0 & adjust_re) {
    
    # Calculate the means of the random effects
    b_means <- apply(b_hat, 2, function(x) do.call(mean_function, c(list(x), mean_function_arguments_random)))
    
    
    # Perform adjustment
    b_hat <- {
      if (k == 1) {
        b_hat - b_means
      } else {
        t(apply(b_hat, 1, "-", b_means))
      } 
    }
    
    beta_hat[X_matching] <- beta_hat[X_matching] + b_means
    
  }
  
  
  rownames(b_hat) <- unique(g)
  
  
  
  results <- list(
    beta = beta_hat,
    beta_init = beta_init,
    b = b_hat,
    sigma = sigma_hat,
    theta = theta_hat
  )
  
  class(results) <- "rankLME"
  
  return(results)
}


res = list()
i=1

for (constant in c(0.7, 0.8, 0.9, 1)) {
  hod <- function(x) hodges_lehmann(x, const = constant)
  Qn_corr2 <- function(x) Qn_corrected(x, p = (4 + n.groups*2))
  
  set.seed(1212)
  r = replicate(50, 
                { 
                  
                  d <- raw_data(
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
                    beta = rep(1, 4)
                  )
                  with(c(d, hod=hod), {
                    ranklme_test(X, y, Z, g, 
                                 intercepts = list(fixed = TRUE, 
                                                   random = TRUE),
                                 mean_function_sd = "hod")
                  }
                  )
                }, simplify = FALSE)
  
  res[[i]] = process_results(r) |> mutate(const = constant)
  i=i+1
}

do.call("rbind", res) |>
  pivot_longer(-c(const, index)) |>
  ggplot(aes(name, value, color=factor(const))) +
  geom_boxplot() +
  geom_hline(yintercept=c(0.5, 1), linetype="dotted")
obj = do.call("rbind", res)

save(obj, file="comparison_regularization_constants.RData")
