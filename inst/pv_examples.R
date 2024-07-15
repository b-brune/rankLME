##############################################################################
#
# Code to reproduce the examples on PV data in Brune, Ortner, Filzmoser (2024)
#
##############################################################################


library(rankLME)
library(tidyverse)
library(lme4)


###
# Example 1: Comparison of Moderate 2 and Moderate 5
###

data(pv1)

mm <- ranklme(
  y = as.matrix(pv1$pmpp_normalized) - 1, 
  X = cbind(pv1$cycle_normalized^2, pv1$cycle_mod5^2, pv1$ramp),
  Z = as.matrix(pv1$cycle_normalized^2), 
  g = as.numeric(as.factor(pv1$name)), 
  adjust_re = TRUE,
  weighted = FALSE, 
  weight_re = FALSE,
  maxit = 150,
  tol = 1e-4,
  intercepts = list(random = FALSE, fixed = FALSE))

mn <- lmer(pmpp_normalized ~ -1 + I(cycle_normalized^2) + I(cycle_mod5^2) + ramp + (0 + I(cycle_normalized^2) | name),
           pv1, offset = rep(1, 48))


cbind(pv1,
      fitted_marginal = cbind(pv1$cycle_normalized^2, pv1$cycle_mod5^2, pv1$ramp) %*% mm$beta + 1,
      fitted_reml = cbind(pv1$cycle_normalized^2, pv1$cycle_mod5^2, pv1$ramp) %*% fixef(mn) + 1) %>%
  rename(module_name = name) %>%
  pivot_longer(c(pmpp_normalized, fitted_marginal, fitted_reml)) %>%
  mutate(name = case_when(name == "pmpp_normalized" ~ "pmpp", 
                          name == "fitted_marginal" ~ "Rank", 
                          name == "fitted_reml" ~ "REML")) %>% 
  mutate(setting = ifelse(setting == "Moderate2", "Moderate 2", "Moderate 5")) %>% 
  ggplot(aes(cycle, value, color = name, shape = name, linetype = name, fill = module_name)) +
  geom_line(alpha = 0.5, linewidth = 1) +
  geom_point(fill="#FA681E") +
  facet_wrap(. ~ setting, scales = "free_x", nrow=1) +
  ylab(expression(paste(P[MPP], " (normalized)"))) + xlab("Treatment cycle") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_color_manual(values = c("pmpp" = "gray40", "REML" = "#458B74", "Rank" = "#FA681E")) +
  scale_linetype_manual(values =  c("pmpp" = "solid", "REML" = "dotdash", "Rank" = "solid")) +
  scale_shape_manual(values = c("pmpp" = 4, "REML" = 16, "Rank" = 25)) 


###
# Example 2: Tropical 1 and Tropical 2 with FTIR peaks
###

data(pv2)

peaks <- c("ft_01area", 
           "ft_03area", 
           "ft_28area")

weighted_rank <- ranklme(
  y = as.matrix(pv2$pmpp_normalized) - 1, 
  X = cbind(as.matrix(pv2 |> dplyr::select(all_of(peaks))), 
            pv2$ramp#, m1_prep$cycle_normalized^2
  ),
  Z = as.matrix(pv2$cycle_normalized^2), 
  g = as.numeric(as.factor(pv2$module_name)), 
  adjust_re = FALSE,
  weighted = TRUE, 
  weight_re = FALSE, 
  mcd = TRUE, 
  maxit = 100,
  control_mean_sd = list(
    mean_function_arguments_fixed = list(),
    mean_function_arguments_random = list(), 
    sd_function_arguments_fixed = list(p = length(peaks) + 1 + length(unique(as.numeric(as.factor(pv2$module_name))))),
    sd_function_arguments_random = list()
  ),
  intercepts = list(random = F, fixed = F),
  leverage_columns = 1:length(peaks),
  tol = 0.0005,
  mean_function = "mean")

weighted_rank

rank <- ranklme(
  y = as.matrix(pv2$pmpp_normalized) - 1, 
  X = cbind(as.matrix(pv2 |> dplyr::select(all_of(peaks))), 
            pv2$ramp#,
            # m1_prep$cycle_normalized^2
  ),
  Z = as.matrix(pv2$cycle_normalized^2), 
  g = as.numeric(as.factor(pv2$module_name)), 
  adjust_re = FALSE,
  weighted = FALSE, 
  weight_re = FALSE, 
  maxit = 100,
  control_mean_sd = list(
    mean_function_arguments_fixed = list(),
    mean_function_arguments_random = list(), 
    sd_function_arguments_fixed = list(p = length(peaks) + 1 + length(unique(as.numeric(as.factor(pv2$module_name))))),
    sd_function_arguments_random = list()
  ),
  intercepts = list(random = F, fixed = F),
  leverage_columns = 1:length(peaks),
  tol = 1e-3,
  mean_function = "mean")


rank

reml <- lme4::lmer(pmpp_normalized ~ -1 + ft_01area +# ft_02area +
                  ft_03area +
                  #ft_08area + ft_10area + ft_23area + 
                  ft_28area + ramp +  #I(cycle_normalized^2) +
                  (0 +  I(cycle_normalized^2) || module_name), pv2, 
                offset = rep(1, nrow(pv2))
)

reml

pv2 %>%
  mutate(weights = weighted_rank$diagnostics$overall_weights) %>%
  rename(`WN 795` = ft_01area,
         #Peak02 = ft_02area,         
         `WN 872` = ft_03area,
         `WN 3426` = ft_28area,
         `pmpp (normalized)` = pmpp_normalized) %>% 
  pivot_longer(c(`pmpp (normalized)`, 
                 `WN 795`, `WN 872`, `WN 3426`), names_to = "peak") %>%
  mutate(peak = factor(peak, levels=c("pmpp (normalized)", "WN 795", "WN 872", "WN 3426"))) %>%
  mutate(Outlying = ifelse(weights < 1, "Yes", "No")) %>%
  ggplot(aes(cycle, value, fill = module_name, linetype=factor(module_number))) +
  geom_line(aes(color = setting)) +
  geom_point(aes(cycle, value, shape = Outlying, color=setting), inherit.aes=FALSE) +
  facet_wrap(peak ~ ., scales = "free_y", nrow=1) +
  scale_color_manual(values = c("Tropical1" = "gray", "Tropical2" = "black")) +
  guides(fill="none", shape=guide_legend(title = "Outlying"), 
         color=guide_legend(title = ""), linetype="none") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_size_identity() +
  scale_linetype_manual(values = c("1" = "solid", "2" = "dotted", "3" = "longdash")) +
  scale_shape_manual(values=c("Yes" = 4, "No" = 16)) +
  xlab("Treatment cycle") + 
  ylab("Measurement") 

