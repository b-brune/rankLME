# diagnostic plots

#' Diagostic plots for rankLME object
#' 
#' @description Provides different diagnostic plots to identify unusual observations
#' 
#' @param x the rankLME object
#' @param which type of plot, where `1` corresponds to ...

plot.rankLME <- function(x, which) {
  
  stopifnot(length(which) == 1)
  
  info = x$diagnostics
  p = length(x$beta) - 1
  
  if (which == 1) {
    ggplot2::ggplot(
      info,
      aes(1:nrow(info), y, color=factor(g))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Index")
  } else if (which == 2) {
    ggplot2::ggplot(
      info,
      aes(fitted_conditional, conditional_residuals, color=factor(g))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = 0, linetype="dotted")
  } else if (which == 3) {
    ggplot2::ggplot(
      info,
      aes(sqrt(overall_leverage), conditional_residuals / x$sigma, color=abs((conditional_residuals / x$sigma)) > 2 | overall_leverage > qchisq(0.95, p))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "none") +
      geom_hline(yintercept = c(-2, 2)) +    
      geom_vline(xintercept = sqrt(qchisq(0.95, p))) +
      xlab("Leverage") + ylab("Standardized conditional residuals")
  } else if (which == 4) {
    
  }
  
}