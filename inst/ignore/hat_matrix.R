library(tidyverse)
library(rankLME)

load("inst/results_leverage_with_hat.RData")
source("inst/simulation_helpers.R")

create_method_column <- function(i) {
  #sapply(x, function(i) {
  if (str_detect(i, "hat_weighted_rank")) {
    "Hat-weighted Rank"
  } else if (str_detect(i, "weighted_rank")) {
    "Weighted Rank"
  } else if (str_detect(i, "rank")) {
    "Rank"
  } else if (str_detect(i, "lme")) {
    "REML"
  } else if (str_detect(i, "koller")) {
    "SMDM"
  } else {
    i
  }
  #}) %>% unname()
}


coefficients_tex <- function(i) {
  #sapply(x, function(i) {
  if (i == "beta_0") {
    "$\\alpha$"
  } else if (i == "beta_1") {
    "$\\beta_1$"
  } else if (i == "beta_2") {
    "$\\beta_2$"
  } else if (i == "beta_3") {
    "$\\beta_3$"
  } else if (i == "sigma") {
    "$\\sigma$"
  } else if (i == "theta_0") {
    "$\\theta_0$"
  } else if (i == "theta_1") {
    "$\\theta_1$"
  } else {
    i
  }
  #}) %>% unname()
}



colorcode <- c("REML" = "#458B74", "Rank" = "#FA681E", "SMDM" = "#0B4768",
               "Weighted Rank" = "#FF69B4", "Hat-weighted Rank" = "gray60")


tikzDevice::tikz("supp_hat_matrix_leverage.tex", width=8, height=4, standAlone=TRUE)
mses |>
  tidyr::pivot_longer(-c(method, outlier_size)) |>
  mutate(method = sapply(method, create_method_column),
         name = sapply(name, function(i) coefficients_tex(i))) |>
  filter(outlier_size %in% c(1, seq(2,20,2))) |>
  ggplot(aes(outlier_size, value, color=method)) +
  geom_line() + geom_point() +
  facet_wrap(name ~ ., scales="free_y", nrow=2) +
  theme_minimal() +
  scale_color_manual(values=colorcode) + theme(legend.position = "bottom") +
  guides(color=guide_legend("")) +
  ylab("MSE") + xlab("Outlier size")
dev.off()

#############

# Simulation with dummy variables:

library(rankLME)


sim <- function(outlier_size) {
  d <- raw_data(
    n = 200, 
    p = 5, 
    k = 2, 
    n.groups = 20, 
    group_size = 10,
    predictors="mixed",
    ranef_distribution = "normal",
    error_distribution = "normal",
    sd_b = 0.5, 
    sd_eps = 1,
    sd_X = 2, 
    beta = rep(1, 5)
  )
  
  d_leverage = leverage(d, outlier_proportion=0.1, outlier_level="single", outlying_columns=c(2,3),
                        single_outlier_type="random", outlier_size=outlier_size, type="variance")
  
  rank = with(d_leverage, ranklme(X, y, Z, g, weighted = FALSE))
  weighted_rank = with(d_leverage, ranklme(X, y, Z, g, weighted = TRUE, weighting_type="malahanobis"))
  hat_weighted_rank = with(d_leverage, ranklme(X, y, Z, g, weighted = TRUE, weighting_type="hat_matrix"))
  
  return(list(
    rank = rank, weighted_rank = weighted_rank, hat_weighted_rank = hat_weighted_rank
  ))
}

set.seed(121212)
res = vector("list", length=length(seq(1,20,2)))
sequence = seq(1, 20, 2)

for (i in seq_along(seq(1, 20, 2))) {
  res[[i]] <- replicate(10, sim(sequence[i]), simplify=FALSE) 
  print(i)
}


get_coefficients <- function(res, coef, method, cols) {
  do.call("rbind", lapply(res, function(x) matrix(unlist(lapply(x, function(y) y[[method]][[coef]])), ncol=cols, byrow=TRUE))) |>
    as.data.frame() |>
    dplyr::mutate(outlier_size = rep(sequence, each=10), method=method)
}


betas = do.call("rbind", lapply(c("rank", "weighted_rank", "hat_weighted_rank"), function(method) get_coefficients(res, "beta", method, 5)) ) 
colnames(betas)[1:5] = c("$\\alpha$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$")

library(tidyverse)

betas |>
  pivot_longer(-c("method", "outlier_size")) |>
  group_by(method, outlier_size, name) |>
  summarize(value = mean((value - 1)^2)) |>
  ggplot(aes(outlier_size, value, color=method)) +
  geom_point() + geom_line() +
  facet_wrap(name ~ .)
  