% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranklme.R
\name{ranklme}
\alias{ranklme}
\title{Rank-based fit for mixed effects models}
\usage{
ranklme(
  X,
  y,
  Z,
  g,
  maxit = 10,
  tol = 1e-06,
  intercepts = list(fixed = TRUE, random = TRUE),
  adjust_re = TRUE,
  weighted = FALSE,
  weight_re = weighted,
  use_outlyingness_weights = TRUE,
  leverage_columns = 1:ncol(X),
  mean_function = "hodges_lehmann",
  sd_function = "Qn_corrected",
  control_mean_sd = list(mean_function_arguments_fixed = list(),
    mean_function_arguments_random = list(), sd_function_arguments_fixed = list(),
    sd_function_arguments_random = list()),
  mcd = TRUE
)
}
\arguments{
\item{X}{(n x p) matrix of regressors, oredered by groups}

\item{y}{(n x 1) target}

\item{Z}{(n x k) matrix with variables corresponding to random effects, usually a subset
of X; subsetting columns should correspond to the first k columns of Z}

\item{g}{(n x 1) vector of group matchings for the observations}

\item{maxit}{maximum number of iterations}

\item{tol}{tolerance until convergence}

\item{intercepts}{a named list of length two with elements \code{fixed}  and \code{random}
indicate whether intercepts should be fitted or not}

\item{adjust_re}{logical, indicates whether the random effects should be adjusted for the mean}

\item{weighted}{logical, should leverage weights be used, defaults to FALSE}

\item{weight_re}{logical, should leverage weights be used for the random slope matrix Z, defaults to `weighted}

\item{use_outlyingness_weights}{logical, usage of outlyingness weights to prevent outliers from contaminating groups}

\item{leverage_columns}{which columns of X should be used for the leverage weights (relevant in cases where we have
random and fixed predictors)}

\item{mean_function}{string, mean function that should be used for centering}

\item{sd_function}{string, scale estimator that should be used for scaling}

\item{control_mean_sd}{list with possible arguments that are handed over to mean and sd function (e.g. cutoff values
for robust estimators)}

\item{mcd}{should the MCD estimator be used to calculate leverage? defaults to TRUE, but can be set to FALSE if
the groups are too small to apply MCD}
}
\value{
A named list with elements
\itemize{
\item \code{beta} -- the estimated fixed effects
\item \code{beta_init} -- the initial value for the fixed effects
\item \code{b} -- the estimated random effects
\item \code{sigma} -- the estimated standard deviations
\item \code{theta} -- the estimated random effects standard deviations
\item \code{diagnostics} -- a data frame with different diagnostic measures, fitted values and residuals
\item some more elements
}
}
\description{
Fits the rank-based mixed effects model as proposed in Brune, Ortner and Filzmoser (2022+).
}
