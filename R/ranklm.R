#' RANK-BASED REGRESSION
#'

#' Fitting the rank based linear model with starting values given by lm
#' Input:
#' @param X regression matrix
#' @param y target
#' @param intercept logical, indicates whether an intercept should be modeled
#' @param weighted logical, should leverage weights be applied
#' @param leverage_columns indicates the columns of X that should be used for calculation
#' of the leverage weights (ignored if weighted = FALSE)
#' @param weighting_type Should malahanobis type weights be used or hat matrix based weights (`malahanobis` or `hat_matrix`)
#' @param mcd should MCD be used for covariance estimation or classic estimates (might be necessary for very small groups)
#' @param ... additional arguments that may be passed to the norm
#'
#' @return The vector of coefficients
#' 
#' @export

ranklm <- function(X, y, intercept = TRUE, 
                   weighted = FALSE,
                   leverage_columns = 1:ncol(X),
                   V = NULL, weights = NULL,
                   c = qchisq(0.95, df = ncol(X)),
                   mean_function = "hodges_lehmann",
                   mean_function_arguments = list(),
                   weighting_type = "malahanobis",
                   mcd = TRUE,
                   ...) {
  
  # If intercept is contained in X, remove it
  
  if (all(X[, 1] == 1)) X <- X[, -1, drop = FALSE]
  
  p <- ncol(X) + ifelse(intercept, 1, 0)
  
  starting_vals <- tryCatch(coef(lm(y ~ X - 1)), error = function(e) NA)
  
  if (any(!is.finite(starting_vals))) starting_vals <- rep(1, ncol(X))
  
  # It is sufficient to calculate V and the weights once:
  
  if (weighted) { 
    
    X_lev <- X[, leverage_columns, drop = FALSE]
    
    if (weighting_type == "malahanobis") {
      
      if (is.null(weights) & !is.null(V)) {
        weights <- sapply(1:nrow(X_lev), function(i) {
          min(1, c / crossprod(X_lev[i, ] - V$center, matpow(V$cov, -1) %*% (X_lev[i, ] - V$center)))
        })
        
      } else if (is.null(V) & is.null(weights)) { 
        
        if (mcd) {
          V <- robustbase::covMcd(X_lev) 
        } else {
          V = list(cov = cov(X_lev), center = colMeans(X_lev))
        }
        
        weights <- sapply(1:nrow(X_lev), function(i) {
          min(1, c / crossprod(X_lev[i, ] - V$center, matpow(V$cov, -1) %*% (X_lev[i, ] - V$center)))
        })
      }
    
    } else if (weighting_type == "hat_matrix") {
      
      hat_matrix <- X_lev %*% solve(t(X_lev) %*% X_lev) %*% t(X_lev)
      
      weights = pmin(1, (2  * ncol(X_lev) / nrow(X_lev)) / diag(hat_matrix))
       
    }
    
  } else {
    
    weights <- rep(1, nrow(X))
    
  }

  
  res <- optim(par = starting_vals, 
               function(beta2) wilcox_norm_weighted(
                 X = X, y = y, beta = beta2, 
                 weighted = weighted, weights = weights, 
                 V = V),
               method = "BFGS")
  
  betahat <- res$par
  
  if (isTRUE(intercept)) {
    
    yhat <- X %*% betahat
    ehat <- y - yhat
    
    # Estimate the intercept
    alphahat <- do.call(mean_function, c(list(ehat), mean_function_arguments)) 
    
    
    yhat <- yhat + alphahat
    ehat <- ehat - alphahat
    
    coef <- c(alphahat, betahat)
    
  } else {
    
    yhat <- X %*% betahat
    
    coef <- betahat
  
  }
  
  final_scores <- a(order(y - yhat))
  
  
    
  return(list(
    coef = coef, 
    final_scores = final_scores, 
    weights = weights)
    )
  
}
