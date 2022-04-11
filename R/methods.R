#' `print method for rankLME
#' 
#' @export

print.rankLME <- function(x, ...) {
  cat("Rank-based mixed effects model", ifelse(x$weighted, "with weighting:", "without weighting:"), "\n\n")
  
  cat("Estimated coefficients:\n")
  cat(unname(x$beta))
  
  cat("\nEstimated variance components:\n\n")
  
  cat("Residual variance:", x$sigma, "\n\n")
  
  cat("Random effects variances:\n")
  cat(x$theta)
  
  cat("\n\nConvergence after", x$iterations, "Iterations.")
}


#' `fitted` method for rankLME
#' 
#' @export

fitted.rankLME <- function(object, type = "conditional", ...) {

  if (type == "conditional") {
    
    return(object$diagnostics$fitted)
  
  } else if (type == "marginal") {

    return(object$diagnostics$y - object$diagnostics$marginal_residuals)
    
  } else {
    
    warning("`type` has to be either \"marginal\" or \"conditional\". Returning conditional fitted values.")
    
    return(object$diagnostics$fitted)
  
  }
  
}


#' `residuals` method for rankLME
#' 
#' @export

residuals.rankLME <- function(object, type = "conditional", ...) {
  
  if (type == "conditional") {
    
    return(object$diagnostics$conditional_residuals)
    
  } else if (type == "marginal") {
    
    return(object$diagnostics$marginal_residuals)
    
  } else {
    
    warning("`type` has to be either \"marginal\" or \"conditional\". Returning conditional residuals.")
    
    return(object$diagnostics$fitted)
    
  }
  
}


