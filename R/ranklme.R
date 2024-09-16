#' Rank-based fit for mixed effects models
#' 
#' @description 
#' Fits the rank-based mixed effects model as proposed in Brune, Ortner and Filzmoser (2022+).
#' 
#' @param X (n x p) matrix of regressors, oredered by groups
#' @param y (n x 1) target
#' @param Z (n x k) matrix with variables corresponding to random effects, usually a subset
#'          of X; subsetting columns should correspond to the first k columns of Z
#' @param g (n x 1) vector of group matchings for the observations 
#' @param maxit maximum number of iterations
#' @param tol tolerance until convergence  
#' @param intercepts a named list of length two with elements `fixed`  and `random`
#'        indicate whether intercepts should be fitted or not
#' @param adjust_re logical, indicates whether the random effects should be adjusted for the mean
#' @param weighted logical, should leverage weights be used, defaults to FALSE
#' @param weight_re logical, should leverage weights be used for the random slope matrix Z, defaults to `weighted`
#' @param weighting_type string, specify the type of leverage weights used, either `malahanobis` or `hat_matrix`.
#' @param use_outlyingness_weights logical, usage of outlyingness weights to prevent outliers from contaminating groups
#' @param leverage_columns which columns of X should be used for the leverage weights (relevant in cases where we have
#' random and fixed predictors)
#' @param mean_function string, mean function that should be used for centering
#' @param sd_function string, scale estimator that should be used for scaling 
#' @param control_mean_sd list with possible arguments that are handed over to mean and sd function (e.g. cutoff values
#' for robust estimators)
#' @param mcd should the MCD estimator be used to calculate leverage? defaults to TRUE, but can be set to FALSE if 
#' the groups are too small to apply MCD
#' 
#' @return 
#' A named list with elements
#' * `beta` -- the estimated fixed effects
#' * `beta_init` -- the initial value for the fixed effects
#' * `b` -- the estimated random effects
#' * `sigma` -- the estimated standard deviations
#' * `theta` -- the estimated random effects standard deviations
#' * `diagnostics` -- a data frame with different diagnostic measures, fitted values and residuals
#' * some more elements 
#' 
#' @export  
#' 
#' 

ranklme <- function(
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
  # Fit the model with random intercept:
  #
  #############################################################################
  
  if (ncol(Z) == 1 & intercepts$random & length(unique(Z[, 1])) == 1) {
    
    if (length(unique(X[, 1])) == 1 & intercepts$fixed) X <- X[, -1, drop = FALSE]
    
    n <- nrow(X)
    p <- ncol(X) + ifelse(intercepts$fixed, 1, 0) 
    k <- 1
    
    Z_matching <- { if (intercepts$fixed) { 1 } else { NULL } }
    
    ### Fit iterative rank estimate
    
    #
    # Fit the first regression
    #
    
    m <- ranklm(X = X, 
                y = y, 
                intercept = intercepts$fixed, 
                weighted = weighted, 
                weighting_type = weighting_type,
                mean_function = mean_function,
                mean_function_arguments = mean_function_arguments_fixed,
                leverage_columns = leverage_columns
    )
    
    # Extract coefficients, predicted values and weights (in case of weighting)
    beta_hat <- m$coef  
    y_hat <- { if (intercepts$fixed) { cbind(1, X) } else { X } } %*% beta_hat
    weights <- m$weights
    
    # Calculate the residuals
    marg_residuals <- y - y_hat
    
    # Predict the random effects by averaging over each group's residuals
    b_hat <- matrix(apply(group_limits, 1, function(x) {
      do.call(mean_function, list(marg_residuals[x[1]:x[2]]))
    }), n.groups, k, byrow = TRUE)
    
    
    # Correct the residuals for the random effect
    cond_residuals <- c(unlist(
      apply(cbind(group_limits, b_hat), 1, function(x) {
        marg_residuals[x[1]:x[2]] - x[3]
      })))
    
    
    # Estimates for variance components
    sigma_hat <- do.call(sd_function, c(list(cond_residuals), 
                                        sd_function_arguments_fixed))
    theta_hat <- do.call(sd_function, c(list(b_hat), 
                                        sd_function_arguments_random))
    
    beta_init <- beta_hat
    
    outlyingness_weights <- { 
      if (use_outlyingness_weights) {
        sapply(2 * sigma_hat / abs(cond_residuals), min, 1)
      } else { rep(1, n) } 
    }
    
    l <- 1
    
    repeat {
      
      # Store "old" set of coefficient estimates
      beta_old <- beta_hat; theta_old <- theta_hat; sigma_old <- sigma_hat; b_hat_old <- b_hat
      
      # Calculate and build the covariance matrix
      Sigma_sq <- matrix(0, n, n)
      for (i in 1:nrow(group_limits)) {
        Sigma_sq[group_limits[i, 1]:group_limits[i, 2], 
                 group_limits[i, 1]:group_limits[i, 2]] <- 
          matpow(
            sigma_hat^2 * diag(group_sizes[i]) +
              theta_hat^2 * matrix(1, group_sizes[i], group_sizes[i]), 
            -1/2)
      }
      
      
      # Rescale y and X with the covariance matrix
      y_star <- Sigma_sq %*% diag(outlyingness_weights) %*% y
      X_star <- Sigma_sq %*% diag(outlyingness_weights) %*% X
      
      # Fit the refined model
      m <- ranklm(
        y = y_star, X = X_star, 
        intercept = intercepts$fixed, 
        weighted = weighted, 
        weighting_type = weighting_type,
        mean_function = mean_function,
        mean_function_arguments = mean_function_arguments_fixed
      )
      
      # new_weights <- m$weights 
      
      beta_hat <- m$coef  
      y_hat_tilde <- X %*% { if (intercepts$fixed) { beta_hat[-1] } else { beta_hat } }
      
      marg_residuals <- y - y_hat_tilde
      
      if (intercepts$fixed) {
        
        alpha <- do.call(mean_function, c(list(marg_residuals), mean_function_arguments_fixed))
        
        marg_residuals <- marg_residuals - alpha
        
        y_hat <- y_hat_tilde + alpha
        
        beta_hat[1] <- alpha
        
      }
      
      b_hat <- matrix(apply(group_limits, 1, function(x) {
        do.call(mean_function, list(marg_residuals[x[1]:x[2]]))
      }), n.groups, k, byrow = T) 
      
      
      cond_residuals <- c(apply(cbind(group_limits, b_hat), 1, function(x) {
        marg_residuals[x[1]:x[2]] - x[3]
      }))
      
      
      sigma_hat <- do.call(sd_function, c(list(cond_residuals), sd_function_arguments_fixed))
      theta_hat <- apply(b_hat, 2, function(x) do.call(sd_function, c(list(x), sd_function_arguments_random)))
      
      outlyingness_weights <- { 
        if (use_outlyingness_weights) {
          sapply(2 * sigma_hat / abs(cond_residuals), min, 1)
        } else { rep(1, n) } }
      
      if (l > maxit | 
          sum((beta_old - beta_hat)^2) / sum(beta_old^2) < tol & 
          sum((c(sigma_hat, theta_hat) - c(sigma_old, theta_old))^2) / sum(c(sigma_old, theta_old)^2) < tol
      ) break
      
      l <- l + 1
      
    }
    
    if (isTRUE(Z_matching) & adjust_re) {
      
      # Calculate the means of the random effects
      b_means <- do.call(mean_function, c(list(b_hat), mean_function_arguments_random))
      
      
      # Perform adjustment
      b_hat <- b_hat - b_means
      beta_hat[Z_matching] <- beta_hat[Z_matching] + b_means
      
    }
    
    ## HIER FEHLT NOCH RECHT VIEL RÜCKGABE STUFF:
    
    rownames(b_hat) <- unique(g)
    
    # fitted <- data.frame(y_hat = y - cond_residuals, g = g)
    # residuals <- data.frame(residuals = cond_residuals, g = g)
    
  } else {
    
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
             mean_function = mean_function,
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
    sigma_hat <- do.call(sd_function, c(list(cond_residuals), sd_function_arguments_fixed))
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
               mean_function = mean_function, mcd = mcd,
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
      
      sigma_hat <- do.call(sd_function, c(list(cond_residuals), sd_function_arguments_fixed))
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
  }
  
  
  ############################################################################
  #
  # Calculate diagnostic measures necessary for plotting:
  #
  ############################################################################
  
  groupwise_variance_diagnostic <- apply(group_limits, 1, function(x) {
    mean((diag(x[2]- x[1] + 1) - Sigma_sq[x[1]:x[2], x[1]:x[2]] %*% 
            tcrossprod(marg_residuals[x[1]:x[2]]) %*% 
            Sigma_sq[x[1]:x[2], x[1]:x[2]])^2)
  })
  
  
  # Calculate leverage of the points:
  X_lev <- X[, leverage_columns, drop = FALSE]
  
  V <- robustbase::covMcd(X_lev)
  lev <- tryCatch(sapply(1:nrow(X_lev), function(i) {
    crossprod(X_lev[i, ] - V$center, matpow(V$cov, -1) %*% (X_lev[i, ] - V$center))
  }), error = function(e) NA)
  
  # Groupwise leverage based on Z-matrix:
  
  lev_groups <- unlist(apply(group_limits, 1, function(x) {
    
    Z_tmp = Z[x[1]:x[2], , drop = FALSE]
    
    if (mcd) {
      V <- robustbase::covMcd(Z_tmp)
    } else {
      V <- list(cov = cov(Z_tmp), center = colMeans(Z_tmp))
    }
    
    sapply(1:nrow(Z_tmp), function(i) {
      crossprod(Z_tmp[i, ] - V$center, matpow(V$cov, -1) %*% 
                  (Z_tmp[i, ] - V$center))
    })
    
  }))
  
  }
  
  
  # Scores in groups: Die sind ein bisschen doof zu berechnen...
  groupwise_scores <-  unlist(lapply(groupwise_models, "[[", 'final_scores'))
  
  if (!weight_re) {
    
    groupwise_models <- apply(group_limits, 1, function(x) {
      ranklm(y = marg_residuals[x[1]:x[2]],
             X = Z[x[1]:x[2], , drop = FALSE],
             intercept = intercepts$random,
             weighted = TRUE, mcd = mcd,
             mean_function = mean_function
      )
    })
    groupwise_leverage_weights <- unlist(lapply(groupwise_models, "[[", 'weights'))
    
  }
  
  # Calculate leverage based on hat matrix
  hii_X = diag(X_lev %*% solve(t(X_lev) %*% X_lev) %*% t(X_lev))
  
  hii_Z <- c(apply(group_limits, 1, function(x) {
    Z_tmp = Z[x[1]:x[2], , drop = FALSE]
    return(diag(Z_tmp %*% solve(t(Z_tmp) %*% Z_tmp) %*% t(Z_tmp)))
  }, simplify=FALSE))
  
  
  # Combine all diagnostic measures into one data frame:
  results_and_diagnostics <- data.frame(
    y = y,
    fitted_conditional = y - cond_residuals, 
    fitted_marginal = y - marg_residuals,
    conditional_residuals = cond_residuals,
    marginal_residuals = marg_residuals,
    g = g,
    var_diagnostic = rep(groupwise_variance_diagnostic, group_sizes),
    overall_leverage = lev,
    overall_hat = hii_X,
    groupwise_leverage = c(lev_groups),
    groupwise_hat = hii_Z,
    overall_scores = a(order(marg_residuals)) ,
    groupwise_scores = groupwise_scores,
    overall_weights = weights, 
    groupwise_weights = unlist(groupwise_leverage_weights),
    outlyingness_weights = sapply(2 * sigma_hat / abs(cond_residuals), min, 1)
  )[order(ord), ]
  
  
  results <- list(
    beta = beta_hat,
    beta_init = beta_init,
    b = b_hat,
    sigma = sigma_hat,
    theta = theta_hat,
    iterations = l,
    diagnostics = results_and_diagnostics,
    weighted = weighted,
    group_limits = group_limits
  )
  
  class(results) <- "rankLME"
  
  return(results)
}