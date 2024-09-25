##############################################################################
#
# Code to reproduce the examples on sleepstudy data in Brune, Ortner, 
#  Filzmoser (2024)
#
##############################################################################

library(lme4)
library(rankLME)
library(robustbase)
library(tidyverse)
library(patchwork)
library(magrittr)

data(sleepstudy)

theme_set(theme_minimal())
colorcode <- c("REML" = "#458B74", "Rank" = "#FA681E", "SMDM" = "#0B4768",
               "Weighted Rank" = "#FF69B4")
ltycode <- c("REML" = "dotdash", "Rank" = "solid", "SMDM" = "dotted", "Weighted Rank" = "dashed")
shapecode <- c("REML" = 16, "Rank" = 25, "SMDM" = 3, "Weighted Rank" = 23) 


#' Contamination of subjects and comparison of model fit
#' Input:
#'  ind ... indices of contaminated observations
#'  size ... size of contamination, defaults to `10`
#'  multi ... logical, multiplicative or additive outliers, defaults to `FALSE`
#' Output:
#' Named list with elements:
#'  `rank` = rank-based model, 
#'  `lme` = Non-robust mixed effects model based on REML, 
#'  `plot` = Plot of the fitted values, 
#'  `data` = the contaminated data set.
contaminated <- function(ind, size = 10, multi = F) {
  
  if (multi) {
    sleepstudy$Reaction[ind] <- sleepstudy$Reaction[ind] * size
  } else  {
    sleepstudy$Reaction[ind] <- sleepstudy$Reaction[ind] + size 
  }
  
  
  y <- as.matrix(sleepstudy$Reaction)
  X <- cbind(1, sleepstudy$Days)
  Z <- cbind(1, sleepstudy$Days)
  g <- sleepstudy$Subject
  
  
  res <- ranklme(X = X, y = y, Z = Z, g = g)
  
  res_lme <- lmer(Reaction ~ Days + (1 + Days || Subject), sleepstudy)
  
  sleepstudy$fitted_lme <- fitted(res_lme)
  sleepstudy$fitted_rank <- X %*% res$beta + c(sapply(unique(sleepstudy$Subject), function(i) {
    tmp <- sleepstudy %>% filter(Subject == i)
    cbind(1, tmp$Days) %*% res$b[as.character(i), ]
  }))
  sleepstudy$outlying <- 1:180 %in% ind
  
  sleepstudy$fitted_lme_marginal <- X %*% fixef(res_lme)
  sleepstudy$fitted_rank_marginal <- X %*% res$beta
  
  
  plt <- ggplot(sleepstudy %>%
                  pivot_longer(c(fitted_lme, fitted_rank, 
                                 fitted_lme_marginal, fitted_rank_marginal)) %>%
                  mutate(marginal = str_detect(name, "marginal"),
                         method = ifelse(
                           marginal, str_extract(name, "_(.*)_"),
                           paste0(str_extract(name, "_(.*)"), "_"))) %>%
                  mutate(method = ifelse(method == "_lme_", "REML", "rank")),
                aes(Days, value, color = method, linetype = marginal, 
                    fill = name)) +
    geom_line(linewidth = 1) +
    geom_point(data = sleepstudy, aes(Days, Reaction), inherit.aes = FALSE) +
    facet_wrap(Subject ~ ., scales = "free") +
    theme(legend.position = "bottom") +
    scale_linetype_manual(values = c("TRUE" = "dotted", "FALSE" = "solid")) +
    guides(linetype = "none", shape = "none") +
    scale_color_brewer(palette = "Set1")
  
  return(list(rank = res, lme = res_lme, plot = plt, data = sleepstudy))
  
}


###
# Models on raw data (without contamination)
###

y <- as.matrix(sleepstudy$Reaction)
X <- cbind(1, sleepstudy$Days)
Z <- cbind(1, sleepstudy$Days)
g <- sleepstudy$Subject

res <- ranklme(X = X, y = y, Z = Z, g = g, leverage_columns=1)

# Rank-based mixed effects model without weighting: 
#   
# Estimated coefficients:
# 252.102 10.62682
# Estimated variance components:
#   
# Residual variance: 16.56628 
# 
# Random effects variances:
#   31.27666 6.554773
# 
# Convergence after 2 Iterations.


res_lme <- lmer(Reaction ~ Days + (1 + Days || Subject), sleepstudy)

# Linear mixed model fit by REML ['lmerMod']
# Formula: Reaction ~ Days + ((1 | Subject) + (0 + Days | Subject))
# Data: sleepstudy
# REML criterion at convergence: 1743.669
# Random effects:
#   Groups    Name        Std.Dev.
# Subject   (Intercept) 25.051  
# Subject.1 Days         5.988  
# Residual              25.565  
# Number of obs: 180, groups:  Subject, 18
# Fixed Effects:
#   (Intercept)         Days  
#        251.41        10.47  


###
# Contamination type 1: A whole group is contaminated
###

set.seed(1007)

a2 <- contaminated(1:10, 3, multi = TRUE)

a2$rank
# Rank-based mixed effects model without weighting: 
#   
#   Estimated coefficients:
#   255.019 10.69506
# Estimated variance components:
#   
#   Residual variance: 17.5316 
# 
# Random effects variances:
#   37.39062 6.554846
# 
# Convergence after 2 Iterations.


a2$lme
#Linear mixed model fit by REML ['lmerMod']
# Formula: Reaction ~ Days + ((1 | Subject) + (0 + Days | Subject))
# Data: sleepstudy
# REML criterion at convergence: 1954.312
# Random effects:
#   Groups    Name        Std.Dev.
# Subject   (Intercept) 119.77  
# Subject.1 Days         14.31  
# Residual               40.66  
# Number of obs: 180, groups:  Subject, 18
# Fixed Effects:
#   (Intercept)         Days  
#        278.54        12.89  


### Variance diagnostic plot:

## Calculate variance diagnostic for REML model:

g <- sleepstudy$Subject

n.groups <- length(unique(g))
group_sizes <- c(table(g))
group_limits <- cbind(cumsum(group_sizes) - group_sizes + 1, cumsum(group_sizes))
Sigma_sq <- matrix(0, nrow(sleepstudy), nrow(sleepstudy))

Z <- cbind(1, sleepstudy$Reaction)

theta_hat <- c(attr(VarCorr(a2$lme)$Subject, "stddev"),
               attr(VarCorr(a2$lme)$Subject.1, "stddev"))
sigma_hat <- attr(VarCorr(a2$lme), "sc")


if (length(theta_hat) > 1) { 
  Sigma_b <- diag(theta_hat^2) 
} else { Sigma_b <- theta_hat^2 * diag(1) }


for (i in 1:nrow(group_limits)) {
  Sigma_sq[group_limits[i, 1]:group_limits[i, 2], 
           group_limits[i, 1]:group_limits[i, 2]] <- 
    rankLME:::matpow(
      sigma_hat^2 * diag(group_sizes[i]) +
        Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE] %*% 
        Sigma_b %*% 
        t(Z[group_limits[i, 1]:group_limits[i, 2], , drop = FALSE]), -1/2)
}

marg_residuals <- a2$data$Reaction - a2$data$fitted_lme_marginal

groupwise_variance_diagnostic_lme <- apply(group_limits, 1, function(x) {
  mean((diag(x[2]- x[1] + 1) - Sigma_sq[x[1]:x[2], x[1]:x[2]] %*% 
          tcrossprod(marg_residuals[x[1]:x[2]]) %*% 
          Sigma_sq[x[1]:x[2], x[1]:x[2]])^2)
})



## Combine results:

vardi <- cbind(
  REML = groupwise_variance_diagnostic_lme, 
  rank = a2$rank$diagnostics$var_diagnostic %>% unique()
) %>%
  as_tibble(rownames = "Individual") 


## Visualize results:

vardi %>%
  rename(Rank = rank) %>%
  pivot_longer(c("REML", "Rank")) %>% 
  ggplot(aes(Individual, value, color = name, shape = name, fill = name)) +
  geom_point() +
  theme_bw() + 
  scale_color_manual(values = colorcode[1:2]) +
  scale_fill_manual(values = colorcode[1:2]) +
  scale_shape_manual(values = shapecode[1:2]) +
  xlab("Group / Individual") + ylab("$d_i$ (log-scale)") +
  theme(axis.text.x = element_text(angle = 90), legend.position="bottom", legend.title = element_blank()) + 
  scale_y_log10()


####
# Contamination type 2: Contamination of a single observation within each subject

a3 <- contaminated((1:nrow(sleepstudy))[seq(5, nrow(sleepstudy), 10)], 3,
                   multi = T)

a3$rank
# Rank-based mixed effects model without weighting: 
#   
#   Estimated coefficients:
#   257.6306 10.49885
# Estimated variance components:
#   
#   Residual variance: 20.91561 
# 
# Random effects variances:
#   33.44222 8.036325
# 
# Convergence after 2 Iterations.

a3$lme
# Linear mixed model fit by REML ['lmerMod']
# Formula: Reaction ~ Days + ((1 | Subject) + (0 + Days | Subject))
# Data: sleepstudy
# REML criterion at convergence: 2371.468
# Random effects:
#   Groups    Name        Std.Dev.
# Subject   (Intercept)   0.0   
# Subject.1 Days          0.0   
# Residual              182.6   
# Number of obs: 180, groups:  Subject, 18
# Fixed Effects:
#   (Intercept)         Days  
# 324.880        6.969  
# optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 


##  Residual plot:

consp <- which(abs(res$diagnostics$conditional_residuals / res$sigma) > 2) 
consp2 <- which(abs(residuals(res_lme)) / attr(VarCorr(res_lme), "sc") > 2) 

a3$data %>%
  mutate(REML = (Reaction - fitted_lme) / attr(VarCorr(a3$lme), "sc"),
         Rank =(Reaction - fitted_rank) / a3$rank$sigma) %>%
  pivot_longer(c("REML", "Rank")) %>%
  mutate(Outlying = ifelse(outlying, "Yes", "No")) %>%
  ggplot(aes(1:length(value)/2, value, color = Outlying, shape=Outlying)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed")+
  facet_wrap(name ~ ., scales="free", ncol=2) +
  ylab("Conditional residuals") + xlab("Index") +
  theme_bw() + 
  theme(legend.position="bottom") 
