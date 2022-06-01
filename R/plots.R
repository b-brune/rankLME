# diagnostic plots

#' Diagostic plots for rankLME object
#' 
#' @description Provides different diagnostic plots to identify unusual observations
#' 
#' @param x the rankLME object
#' @param which type of plot, see Details
#' 
#' @details
#' The function creates different diagnostic plots for the rankLME model. The type of 
#' plot is specified by the `which` argument:
#' 
#' * `1` Plot of index vs. y by group color
#' * `2` Plot of fitted values vs. conditional residuals
#' * `3` Plot of overall leverage vs. standardized conditional residuals
#' * `4` Plot of groupwise leverage vs. standardized conditional residuals, faceted by group
#' * `5` Boxplots of the marginal residuals by group
#' * `6` Variance diagnostic
#' 
#' @export

plot.rankLME <- function(x, which=2) {
  
  stopifnot(length(which) == 1)
  
  info = x$diagnostics
  p = length(x$beta) - 1
  k = ncol(x$b) - 1
  
  ggplot2::theme_set(ggplot2::theme_bw())
  
  if (which == 1) {
    ggplot2::ggplot(
      info,
      ggplot2::aes(1:nrow(info), y, color=factor(g))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Index")
  } else if (which == 2) {
    ggplot2::ggplot(
      info,
      ggplot2::aes(fitted_conditional, conditional_residuals, color=factor(g))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = 0, linetype="dotted")
  } else if (which == 3) {
    ggplot2::ggplot(
      info,
      ggplot2::aes(sqrt(overall_leverage), conditional_residuals / x$sigma, color=abs((conditional_residuals / x$sigma)) > 2 | overall_leverage > qchisq(0.95, p))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = c(-2, 2)) +    
      ggplot2::geom_vline(xintercept = sqrt(qchisq(0.95, p))) +
      ggplot2::xlab("Leverage") + 
      ggplot2::ylab("Standardized conditional residuals")
  } else if (which == 4) {
    ggplot2::ggplot(
      info,
      ggplot2::aes(sqrt(groupwise_leverage), conditional_residuals / x$sigma, color=abs((conditional_residuals / x$sigma)) > 2 | groupwise_leverage > qchisq(0.95, k))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype="dotted") +    
      ggplot2::geom_vline(xintercept = sqrt(qchisq(0.95, k)), linetype="dotted") +
      ggplot2::xlab("Group-specific leverage") + 
      ggplot2::ylab("Standardized conditional residuals") +
      ggplot2::facet_wrap(g ~ ., labeller = "label_both")
  } else if (which == 5) {
    ggplot2::ggplot(
      info,
      ggplot2::aes(factor(g), marginal_residuals, color=factor(g))
    ) +
      ggplot2::geom_boxplot() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Group") + ggplot2::ylab("Marginal residuals") +
      ggplot2::geom_hline(yintercept = 0, linetype="dotted") 
  } else if (which == 6) {
    ggplot2::ggplot(
      info,
      ggplot2::aes(factor(g), var_diagnostic, color=factor(g))
    ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position="none") +
      ggplot2::xlab("Group") + ggplot2::ylab("Variance diagnostic measure") 
  }
  
}
