##############################################################################
#
# Code to reproduce the examples on PV data in Brune, Ortner, Filzmoser (2024)
#
##############################################################################

###
# Example 1: Comparison of Moderate 2 and Moderate 5
###

library(rankLME)
library(ggplot)

data(pv1)

mm <- ranklme(
  y = as.matrix(m1_prep$pmpp_normalized) - 1, 
  X = cbind(m1_prep$cycle_normalized^2, m1_prep$cycle_mod5^2, m1_prep$ramp),
  Z = as.matrix(m1_prep$cycle_normalized^2), 
  g = as.numeric(as.factor(m1_prep$name)), 
  adjust_re = TRUE,
  weighted = FALSE, 
  weight_re = FALSE,
  maxit = 150,
  tol = 1e-4,
  intercepts = list(random = FALSE, fixed = FALSE))

mn <- lmer(pmpp_normalized ~ -1 + I(cycle_normalized^2) + I(cycle_mod5^2) + ramp + (0 + I(cycle_normalized^2) | name),
           m1_prep, offset = rep(1, 48))


cbind(m1_prep,
      fitted_marginal = cbind(m1_prep$cycle_normalized^2, m1_prep$cycle_mod5^2, m1_prep$ramp) %*% mm$beta + 1,
      fitted_reml = cbind(m1_prep$cycle_normalized^2, m1_prep$cycle_mod5^2, m1_prep$ramp) %*% fixef(mn) + 1) %>%
  rename(module_name = name) %>%
  pivot_longer(c(pmpp_normalized, fitted_marginal, fitted_reml)) %>%
  mutate(name = case_when(name == "pmpp_normalized" ~ "pmpp", 
                          name == "fitted_marginal" ~ "Rank", 
                          name == "fitted_reml" ~ "REML")) %>% 
  mutate(setting = ifelse(setting == "Moderate2", "Moderate 2", "Moderate 5")) %>% 
  ggplot(aes(cycle, value, color = name, shape = name, linetype = name, fill = module_name)) +
  geom_line(alpha = 0.5, size = 1) +
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

