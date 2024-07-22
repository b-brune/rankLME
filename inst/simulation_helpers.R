# Simulation helper functions

#'
#' Takes a data object (output from `raw_data()` or one of the contamination functions)
#' and fits the models for the comparison (REML, SMDM, rank, weighted_rank)
#' Input:
#' @param d ... The dataset 
#' @param with_smdm ... logical, should the SMDM model be fit as well? (Computationally expensive)
#' 
#' Output:
#' A named list with the fitted models.
#' 
#' 
compare_model_fits <- function(d, with_smdm=FALSE) {
  require(lme4)
  require(robustlmm)
  require(rankLME)
  
  p <- ncol(d$X)
  k <- ncol(d$Z)
  n <- nrow(d$X)
  
  # Make a data matrix
  dat <- with(
    d,
    data.frame(
      magrittr::set_colnames(
        cbind(y, X, Z, g),
        c(
          "y", paste0("X", 1:p - 1),
          paste0("Z", 1:k - 1), "g"
        )
      )
    )
  )
  
  # Build formulas for the functions that need formulas:
  f1 <- paste("y", "~", paste(colnames(dat)[3:(p + 1)], collapse = " + "))
  f2 <- paste("(", "1 +", paste(colnames(dat)[(p + 3):(p + k + 1)], collapse = "+"), "|| g )", collapse = "")
  
  f <- formula(paste(f1, f2, sep = " + "))
  
  weighted_rank <- ranklme(
    X = d$X, y = d$y, Z = d$Z, g = d$g, 
    sd_function = "Qn_corrected",
    adjust_re =TRUE,
    control_mean_sd = list(
      mean_function_arguments_fixed = list(),
      mean_function_arguments_random = list(), 
      sd_function_arguments_fixed = list(p = (p + d$n.groups)),
      sd_function_arguments_random = list()),
    weighted = TRUE, 
    maxit = 20
  )
  
  rank <- ranklme(
    X = d$X, y = d$y, Z = d$Z, g = d$g,
    sd_function = "Qn_corrected", adjust_re =TRUE,
    control_mean_sd = list(
      mean_function_arguments_fixed = list(),
      mean_function_arguments_random = list(), 
      sd_function_arguments_fixed = list(p = (p + d$n.groups)),
      sd_function_arguments_random = list()),    weighted = FALSE,
    maxit = 20
  )
  
  lme <- lmer(f, data = dat)
  
  result_list = list(
    rank = rank, weighted_rank = weighted_rank, lme = lme
  )
  
  if (with_smdm) {
    koller <- rlmer(f, data = dat)
    result_list[[4]] <- koller
    names(result_list)[4] <- "koller"
  }
  
  return(
    result_list
  )
}


#'
#' Takes the fitted models (output from `compare_model_fits()` and extracts the coefficients).
#' Input:
#' @param fitted_object ... a list of models
#' 
#' Output:
#' A data frame with the algorithms as columns, and names of the coefficients as rows.
#' 
#' 
get_coefficients_from_model <- function(fitted_object) {
  lme = with(fitted_object, 
             list(beta = fixef(lme),
                  variances = c(sigma(lme), 
                                unlist(lapply(VarCorr(lme), function(z) attr(z, "stddev")))))
  )
  
  if ("koller" %in% names(fitted_object)) {
    smdm = with(fitted_object,
                list(beta = fixef(koller),
                     variances = c(sigma(koller), unlist(lapply(VarCorr(koller), function(z) attr(z, "stddev")))))
    )
  } else {
    SMDM = NULL
  }
  
  rank = with(fitted_object,
              list(beta = rank$beta,
                   variances = c(rank$sigma, rank$theta))
  )
  
  weighted_rank = with(fitted_object,
                       list(beta = weighted_rank$beta,
                            variances = c(weighted_rank$sigma, weighted_rank$theta))
  )
  
  return(
    data.frame(
      REML = unlist(lme), 
      SMDM = unlist(smdm), 
      rank = unlist(rank), 
      weighted_rank = unlist(weighted_rank), 
      coef = c("beta_0", "beta_1", "beta_2", "beta_3", "error_variance", "randef_variance1", "randef_variance2"))  
  )
}


#'
#' Takes the fitted models (output from `compare_model_fits()` and extracts the coefficients).
#' Input:
#' @param fitted_object ... a list of models
#' 
#' Output:
#' A data frame with the algorithms as rows, and names of the coefficients as columns.
#' 
#' 
get_coefficients_from_simulation <- function(fitted_object) {
  lme = unlist(with(fitted_object, 
                    list(beta = fixef(lme),
                         variances = c(sigma(lme), 
                                       unlist(lapply(VarCorr(lme), function(z) attr(z, "stddev")))))
  ))
  
  if ("koller" %in% names(fitted_object)) {
    smdm = with(fitted_object,
                list(beta = fixef(koller),
                     variances = c(sigma(koller), unlist(lapply(VarCorr(koller), function(z) attr(z, "stddev")))))
    )
  } else {
    smdm = NULL
  }
  
  rank = unlist(with(fitted_object,
                     list(beta = rank$beta,
                          variances = c(rank$sigma, rank$theta))
  ))
  
  weighted_rank = unlist(with(fitted_object,
                              list(beta = weighted_rank$beta,
                                   variances = c(weighted_rank$sigma, weighted_rank$theta))
  ))
  
  
  names(lme) <- names(rank) <- names(weighted_rank) <- c("beta_0", "beta_1", "beta_2", "beta_3", "sigma^2", "theta_0", "theta_1")
  if ("koller" %in% names(fitted_object)) {
    names(smdm) <- names(lme)
  }
  
  if (!("koller" %in% names(fitted_object))) {
    
    result =   dplyr::bind_rows(
      REML = unlist(lme), 
      rank = unlist(rank), 
      weighted_rank = unlist(weighted_rank) 
    ) |>
      dplyr::mutate(method = c(
        "REML", 
        "rank", 
        "weighted_rank"))
  } else {
    result =   dplyr::bind_rows(
      REML = unlist(lme), 
      SMDM = unlist(smdm), 
      rank = unlist(rank), 
      weighted_rank = unlist(weighted_rank) 
    ) |>
      dplyr::mutate(method = c(
        "REML", 
        "SMDM",
        "rank", 
        "weighted_rank"))
  }
  
  return(
    result
  )
  
}



#'
#' Carries out one simulation run and returns the results as a dataframe.
#' The data is contaminated with leverage outliers. The parameters within the
#' function can be adjusted to test out other settings.
#' 
#' Input:
#' @param outlier_size ... The desired outlier size
#' @param with_smdm ... logical, should the SMDM model be fit as well? (Computationally expensive)
#' 
#' Output:
#' A data frame with the algorithms as rows, and names of the coefficients as columns. (output of
#' `get_coefficients_from_simulation()` function with an extra column for outlier_size.)
#' 
replicate_model_fits = function(outlier_size, with_smdm=FALSE) {
  
  d <- raw_data(
    n = 200, 
    p = 4, 
    k = 2, 
    n.groups = 20, 
    group_size = 10,
    ranef_distribution = "normal",
    error_distribution = "normal",
    sd_b = 0.5, 
    sd_eps = 1,
    sd_X = 2, 
    beta = rep(1, 4)
  )
  
  d_contaminated <- leverage(
    d, 
    outlier_proportion=0.1, 
    outlier_level="single", 
    single_outlier_type="random", 
    outlier_size=outlier_size, 
    type="variance"
  )
  
  fitted = compare_model_fits(d_contaminated, with_smdm=with_smdm)
  coefficients = get_coefficients_from_simulation(fitted) |> dplyr::mutate(outlier_size = outlier_size)
  
  return(coefficients)
}
