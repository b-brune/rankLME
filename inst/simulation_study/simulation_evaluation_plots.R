library(glue)

# Set path here:
PATH_TO_SIMULATION_RESULTS = "../../../simulation_results_rankLME"

# Set wd here:
setwd("inst/simulation_study")

source("setup_simulation_evaluation.R")


# ---------------------------------------------
# Bias and efficiency
# ---------------------------------------------

load(glue("{PATH_TO_SIMULATION_RESULTS}/consistency.Rdata"))

consistency = preprocess_simulation_results(consistency)

gc()

###
###
# Boxplots with consistency
###
###

boxplots_consistency(consistency, method="Rank", additional_filters = "sd_function == 'Qn_corrected'") 
boxplots_consistency(consistency, method="Weighted Rank", additional_filters = "sd_function == 'Qn_corrected'") 
boxplots_consistency(consistency, method="SMDM", additional_filters = "sd_function == 'Qn_corrected'") 
boxplots_consistency(consistency, method="REML", additional_filters = "sd_function == 'Qn_corrected'") 


#########
######### CALCULATE AREs
#########


rank <- consistency %>%
  filter(method == "Rank", sd_function=="Qn_corrected") |>
  group_by(n.groups, group_size) %>% 
  select(starts_with("beta") & !contains("init")) %>%
  nest(data = contains("beta")) %>%
  mutate(cov_rank = lapply(data, \(x) cov(as.matrix(x)))) |>
  select(-data)

init <- consistency %>%
  filter(method == "Rank", sd_function=="Qn_corrected") |>
  group_by(n.groups, group_size) %>% 
  select(starts_with("beta") & ends_with("init")) %>%
  nest(data = contains("beta")) %>%
  mutate(cov_init = lapply(data, \(x) cov(as.matrix(x)))) |>
  select(-data)

weighted_rank <- consistency %>%
  filter(method == "Weighted Rank", sd_function=="Qn_corrected") |>
  group_by(n.groups, group_size) %>% 
  select(starts_with("beta") & !ends_with("init")) %>%
  nest(data = contains("beta")) %>%
  mutate(cov_weighted = lapply(data, \(x) cov(as.matrix(x))))|>
  select(-data)


lme <- consistency %>%
  filter(method == "REML") |>
  group_by(n.groups, group_size) %>% 
  select(starts_with("beta") & !ends_with("init")) %>%
  nest(data = contains("beta")) %>%
  mutate(cov_lme = lapply(data, \(x) cov(as.matrix(x))))|>
  select(-data)

koller <- consistency %>%
  filter(method == "SMDM") |>
  group_by(n.groups, group_size) %>% 
  select(starts_with("beta") & !ends_with("init")) %>%
  nest(data = contains("beta")) %>%
  mutate(cov_koller = lapply(data, \(x) robustbase::covMcd(as.matrix(x))$cov))|>
  select(-data)

are <- rank |>
  left_join(init, by = c("n.groups", "group_size")) |>
  ungroup() |>
  mutate(rank_vs_init = sapply(1:nrow(rank), \(i)  (det(cov_rank[[i]]) / det(cov_init[[i]]))^(1/4))) |>
  left_join(lme, by = c("n.groups", "group_size")) |>
  mutate(rank_vs_lme = sapply(1:nrow(rank), \(i)  (det(cov_rank[[i]]) / det(cov_lme[[i]]))^(1/4))) |>
  left_join(weighted_rank, by = c("n.groups", "group_size")) |>
  left_join(koller, by = c("n.groups", "group_size")) |>
  mutate(rank_vs_koller = sapply(1:nrow(rank), \(i)  (det(cov_rank[[i]]) / det(cov_koller[[i]]))^(1/4))) |>
  mutate(rank_vs_weighted = sapply(1:nrow(rank), \(i)  (det(cov_rank[[i]]) / det(cov_weighted[[i]]))^(1/4))) |>
  select(-starts_with("cov")) |>
  set_colnames(c(
    "$g$", "$n_i$", "initial value", "REML", "SMDM", "Weighted Rank"
  ))

# Table with efficiencies:
xtable::xtable(are |> filter(`$n_i$` %in% c(5, 20, 50), `$g$` %in% c(5, 20, 50)), digits = c(0,0,0,3,3,3,3))




# ---------------------------------------------
# Response outliers
# ---------------------------------------------

load(glue("{PATH_TO_SIMULATION_RESULTS}/youtliers.Rdata"))
head(youtliers)

youtliers = preprocess_simulation_results(youtliers)

gc()


# Reaction to multiplicative y-outliers of increasing size (random)
(youtliers %>% filter(setting == 1, n == 400, outlier_proportion == 0.1) %>%
    lineplot_mean(
      includelme = TRUE, 
      additional_filters = "y_outlier_type == 'variance' & outlier_size <= 300",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
      scales = "fixed", variables = "coefficients")) +
  (youtliers %>% filter(setting == 1, n == 400, outlier_proportion == 0.1) %>%
     lineplot_mean(
       includelme = TRUE, 
       additional_filters = "y_outlier_type == 'variance' & outlier_size <= 300",
       aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
       scales = "fixed", variables = "variances")) +
  plot_layout(guides = "collect", nrow=2) & 
  theme(title = element_text(size = 8), axis.text.x = element_text(angle = 90)) & 
  scale_y_log10(labels = scales::comma) &
  coord_cartesian(ylim = c(NA, 30))


# Reaction to additive y-outliers of increasing size (random)
(youtliers %>% filter(setting == 1, n == 400, outlier_proportion == 0.1) %>%
    lineplot_mean(
      includelme = TRUE, 
      additional_filters = "y_outlier_type == 'mean' & outlier_size <= 300",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
      scales = "fixed", variables = "coefficients")) +
  (youtliers %>% filter(setting == 1, n == 400, outlier_proportion == 0.1) %>%
     lineplot_mean(
       includelme = TRUE, 
       additional_filters = "y_outlier_type == 'mean' & outlier_size <= 300",
       aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
       scales = "fixed", variables = "variances")) +
  plot_layout(guides = "collect", nrow=2) & 
  theme(title = element_text(size = 8), axis.text.x = element_text(angle = 90)) & 
  scale_y_log10(labels = scales::comma) &
  coord_cartesian(ylim = c(NA, 30))


# Reaction to multiplicative y-outliers of increasing size: (sequential)
(youtliers %>% 
    filter(setting == 1, n == 400, outlier_proportion == 0.1) %>%
    lineplot_mean(
      includelme = TRUE,
      outlier_positioning = "sequential",
      additional_filters = "y_outlier_type == 'variance' & outlier_size <= 300",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
      scales = "fixed", variables = "coefficients")) +
  (youtliers %>% 
     filter(setting == 1, n == 400, outlier_proportion == 0.1, single_outlier_type == "sequential") %>%
     lineplot_mean(
       includelme = TRUE, 
       outlier_positioning = "sequential",
       additional_filters = "y_outlier_type == 'variance' & outlier_size <= 300",
       aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
       scales = "fixed", variables = "variances")) +
  plot_layout(guides = "collect", nrow=2) & 
  theme(title = element_text(size = 8), axis.text.x = element_text(angle = 90)) & 
  scale_y_log10(labels = scales::comma) &
  coord_cartesian(ylim = c(NA, 30))


# Reaction to additive y-outliers of increasing size: (sequential)
(youtliers %>% filter(setting == 1, n == 400, outlier_proportion == 0.1) %>%
    lineplot_mean(
      includelme = TRUE, 
      outlier_positioning = "sequential",
      additional_filters = "y_outlier_type == 'mean' & outlier_size <= 300",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
      scales = "fixed", variables = "coefficients")) +
  (youtliers %>% filter(setting == 1, n == 400, outlier_proportion == 0.1) %>%
     lineplot_mean(
       includelme = TRUE,
       outlier_positioning = "sequential",
       additional_filters = "y_outlier_type == 'mean' & outlier_size <= 300",
       aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE (log-scale)",
       scales = "fixed", variables = "variances")) +
  plot_layout(guides = "collect", nrow=2) & 
  theme(title = element_text(size = 8), axis.text.x = element_text(angle = 90)) & 
  scale_y_log10(labels = scales::comma) &
  coord_cartesian(ylim = c(NA, 30))


# Efficiency:
(
  plot_mse_efficiency(
    dat = youtliers |> 
      filter(outlier_proportion==0.05, outlier_size<=300, 
             y_outlier_type=="variance", single_outlier_type == "random"), "Outlier proportion: 0.05") +
    plot_mse_efficiency(
      dat = youtliers |> 
        filter(outlier_proportion==0.1, outlier_size<=300, 
               y_outlier_type=="variance", single_outlier_type == "random"), "Outlier proportion: 0.1") ) +
  plot_layout(guides="collect", nrow=2)



# ---------------------------------------------
# "Breakdown point"
# ---------------------------------------------

load(glue("{PATH_TO_SIMULATION_RESULTS}/breakdown.Rdata"))
head(breakdown)

breakdown = preprocess_simulation_results(breakdown)
gc() 


# Breakdown point: Increasing proportion of outliers of size 100
(
  breakdown %>% filter(n == 400) %>%
    filter(outlier_proportion %in% seq(0, 0.5, 0.05)) %>%
    lineplot_mean(includelme = TRUE, 
                  aggregation_function = "calc_mse",
                  xlab = "Outlier proportion",
                  xaxis = "outlier_proportion", 
                  additional_filters = "y_outlier_type == 'variance'",
                  add_hline = FALSE, ylab = "MSE (log-scale)", variables = "coefficients", scales = "fixed")  +
    scale_y_log10()) +
  (breakdown %>% filter(n == 400)  %>% 
     filter(outlier_proportion %in% seq(0, 0.5, 0.05)) %>%
     lineplot_mean(includelme = TRUE, 
                   aggregation_function = "calc_mse", 
                   xlab = "Outlier proportion",
                   xaxis = "outlier_proportion", 
                   additional_filters = "y_outlier_type == 'variance'",
                   add_hline = FALSE, ylab = "MSE (log-scale)", variables = "variances", scales = "fixed")  +
     scale_y_log10()
  ) + 
  plot_layout(guides = "collect", nrow = 2) 



# ---------------------------------------------
# "Breakdown point" (leverage)
# ---------------------------------------------

load(glue("{PATH_TO_SIMULATION_RESULTS}/breakdown_leverage.Rdata"))
head(breakdown_leverage)

breakdown_leverage = preprocess_simulation_results(breakdown_leverage)
gc() 



(breakdown_leverage %>% 
    filter(n == 400, single_outlier_type=="random", outlier_size==100, outlier_proportion <= 0.4) %>%
    lineplot_mean(includelme = TRUE, 
                  aggregation_function = "calc_mse",
                  xlab = "Outlier proportion",
                  xaxis = "outlier_proportion", 
                  additional_filters = "leverage_type == 'variance'",
                  add_hline = FALSE, ylab = "MSE", variables = "coefficients", scales = "free")) +
  (breakdown_leverage %>% filter(n == 400, leverage_type=="variance", outlier_size==100, outlier_proportion <= 0.4)  %>% 
     lineplot_mean(includelme = TRUE, 
                   aggregation_function = "calc_mse", 
                   xlab = "Outlier proportion",
                   xaxis = "outlier_proportion", 
                   additional_filters = "leverage_type == 'variance'",
                   add_hline = FALSE, ylab = "MSE", variables = "variances", scales = "free") 
   ) + 
  plot_layout(guides = "collect", nrow = 2) 



# ---------------------------------------------
# Leverage points
# ---------------------------------------------

load(glue("{PATH_TO_SIMULATION_RESULTS}/leverage.Rdata"))
leverage = preprocess_simulation_results(leverage)
head(leverage)

gc() 


# Reaction to multiplicative leverage points of increasing size (random)
(leverage %>% filter(n == 400, outlier_proportion == 0.1, single_outlier_type == "random") %>%
    lineplot_mean(
      includelme = TRUE, 
      additional_filters = "leverage_type == 'variance'",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
      # title = "Outlier proportion: 0.10", 
      scales = "free", variables = "coefficients")) +
(leverage %>% filter(n == 400, outlier_proportion == 0.1, single_outlier_type == "random") %>%
    lineplot_mean(
      includelme = TRUE, 
      additional_filters = "leverage_type == 'variance'",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
      # title = "Outlier proportion: 0.10", 
      scales = "free", variables = "variances")) +
  plot_layout(guides = "collect", nrow = 2) 


# Reaction to multiplicative leverage points of increasing size (sequential)
(leverage %>% filter(n == 400, outlier_proportion == 0.1) %>%
    lineplot_mean(
      includelme = TRUE, 
      outlier_positioning = "sequential",
      additional_filters = "leverage_type == 'variance'",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
      scales = "free", variables = "coefficients")) +
  (leverage %>% filter(n == 400, outlier_proportion == 0.1) %>%
     lineplot_mean(
       includelme = TRUE, 
       outlier_positioning = "sequential",
       additional_filters = "leverage_type == 'variance'",
       aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
       scales = "free", variables = "variances")) +
  plot_layout(guides = "collect", nrow = 2) 


# Reaction to additive leverage points of increasing size:
(leverage %>% filter(n == 400, outlier_proportion == 0.1) %>%
    lineplot_mean(
      includelme = TRUE, 
      additional_filters = "leverage_type == 'mean'",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
      # title = "Outlier proportion: 0.10", 
      scales = "free", variables = "coefficients")) +
  (leverage %>% filter(n == 400, outlier_proportion == 0.1) %>%
     lineplot_mean(
       includelme = TRUE, 
       additional_filters = "leverage_type == 'mean'",
       aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
       # title = "Outlier proportion: 0.10", 
       scales = "free", variables = "variances")) +
  plot_layout(guides = "collect", nrow = 2) 


# Reaction to additive leverage points of increasing size (sequential)
(leverage %>% filter(n == 400, outlier_proportion == 0.1) %>%
    lineplot_mean(
      includelme = TRUE, 
      outlier_positioning = "sequential",
      additional_filters = "leverage_type == 'mean'",
      aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
      # title = "Outlier proportion: 0.10", 
      scales = "free", variables = "coefficients")) +
  (leverage %>% filter(n == 400, outlier_proportion == 0.1) %>%
     lineplot_mean(
       includelme = TRUE, 
       outlier_positioning = "sequential",
       additional_filters = "leverage_type == 'mean'",
       aggregation_function = "calc_mse", add_hline = FALSE, ylab = "MSE",
       # title = "Outlier proportion: 0.10", 
       scales = "free", variables = "variances")) +
  plot_layout(guides = "collect", nrow = 2) 






