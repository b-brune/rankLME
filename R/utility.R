##
##
## UTILITY FUNCTIONS
##
##


#' Generalized matrix power function
#' Implemented by Daniel Kapla
#' (borrowed from my tvRRR package)
#' 

matpow <- function (A, pow, symmetric = FALSE, tol = 1e-07) {
  if (nrow(A) != ncol(A)) {
    stop("Expected a square matix, but 'A' is ", nrow(A), 
         " by ", ncol(A))
  }
  if (pow > 0) {
    if (pow == 1) {
      return(A)
    }
    svdA <- La.svd(A)
    return(svdA$u %*% ((svdA$d^pow) * svdA$vt))
  }
  else if (pow == 0) {
    return(diag(nrow(A)))
  }
  else {
    qrA <- qr(A, tol = tol)
    if (qrA$rank == nrow(A)) {
      return(matpow(solve(qrA), abs(pow), tol = tol))
    }
    else {
      svdA <- svd(A)
      positives <- svdA$d > tol * svdA$d[1]
      d <- c(svdA$d[positives]^pow, rep(0, sum(!positives)))
      if (symmetric) {
        return(svdA$v %*% (d * t(svdA$v)))
      }
      else {
        return(svdA$v %*% (d * t(svdA$u)))
      }
    }
  }
}



#' Wilcoxon scores 
#' 
#' Input: 
#' @param v a numeric vector
#' 
#' Output:
#' A numeric vector of the Wilcoxon scores

a <- function(v) {
  
  u <- order(v)
  u <- u / (length(v) + 1)
  a <- sqrt(12) * (u - 1/2)
  
  return(a)
}

#' Wilcox norm function
wilcox_norm <- function(X, y, beta) {
  sum(a(order(y - X %*% beta)) * (y - X %*% beta))
}

#' Weighted Wilcoxon norm, following 2004 Paper by McKean
#'
#' @examples
#' X <- replicate(4, rnorm(100))
#' beta <- c(1, 2, 3, 4)
#' y <- X %*% beta + rnorm(100)
#' X[1:10, ] <- 50
#' wilcox_norm_weighted(X, y, beta, weighted = TRUE)

wilcox_norm_weighted <- function(X, y, beta, 
                                 weighted = TRUE, 
                                 c = qchisq(0.95, df = ncol(X)),
                                 weights = NULL,
                                 V = NULL) {

  if (!weighted) {
  
    return( sum(a(order(y - X %*% beta)) * (y - X %*% beta)) )
    
  } else if (weighted) {
    
    if (is.null(weights)) { 
      
      if (is.null(V)) { V <- robustbase::covMcd(X) }
          
      weights <- sapply(1:nrow(X), function(i) {
        min(1, c / crossprod(X[i, ] - V$center, matpow(V$cov, -1) %*% (X[i, ] - V$center)))
      })
    } 
    
    w_mat <- outer(weights, weights)
    
    return(sum(w_mat[lower.tri(w_mat)] * as.vector(dist(y - X %*% beta, method = "manhattan"))))
    
  }
}
  
# X <- replicate(4, rnorm(100))
# beta <- c(1, 2, 3, 4)
# y <- X %*% beta + rnorm(100)
# X[1:10, ] <- 50
# 
# wilcox_norm_weighted(X, y, beta, weighted = TRUE)


#' HL1 estimator (from robnptests)
hodges_lehmann <- function (x, na.rm = FALSE) {
  
  if (!na.rm & any(is.na(x))) {
    return(NA_real_)
  }
  else if (na.rm & any(is.na(x))) {
    x <- as.vector(stats::na.omit(x))
  }
  x.grid <- cbind(rep(seq_along(x), each = length(x)), seq_along(x))
  x.diffs <- x.grid[x.grid[, 1] < x.grid[, 2], , drop = FALSE]
  mean.pairwise.sums <- (x[x.diffs[, 1]] + x[x.diffs[, 2]])/2
  return(stats::median(mean.pairwise.sums))
}



#' RLME dispvar
dispvar <- function (x, p = 0) {
  
  n = length(x)
  rx = rank(x, ties.method = c("random"))
  
  sc = sqrt(12) * ((rx/(n + 1)) - 0.5)
  dispvar = sqrt(pi/3) * sum(x * sc) / (n - p)
  
  dispvar
}



#' Qn-variance corrected

Qn_corrected <- function(x, p = 0, ...) {
  n <- length(x)
  sqrt(n / (n - p)) * robustbase::Qn(x, ...)
} 

#' Univariate MCD estimator
uniMCD <- function(x, p = 0, ...) {
  sqrt(n / (n - p) * robustbase::covMcd(x)$cov)
} 
