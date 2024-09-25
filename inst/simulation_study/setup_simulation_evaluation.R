# Functions to process simulation results
library(tidyverse)
library(data.table)
library(magrittr)
library(gridExtra)
library(patchwork)
library(glue)

theme_custom <- function() { 
  #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      legend.position = "bottom"
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}

coefficients_tex <- function(x) {
  sapply(x, function(i) {
    if (str_detect(i, "beta1") | str_detect(i, "beta_init1")) {
      "$\\alpha$"
    } else if (str_detect(i, "beta2") | str_detect(i, "beta_init2")) {
      "$\\beta_1$"
    } else if (str_detect(i, "beta3") | str_detect(i, "beta_init3")) {
      "$\\beta_2$"
    } else if (str_detect(i, "beta4") | str_detect(i, "beta_init4")) {
      "$\\beta_3$"
    } else if (str_detect(i, "sigma")) {
      "$\\sigma$"
    } else if (str_detect(i, "theta0")) {
      "$\\theta_0$"
    } else if (str_detect(i, "theta1")) {
      "$\\theta_1$"
    }
  }) %>% unname()
}


create_method_column <- function(x) {
  sapply(x, function(i) {
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
    }
  }) %>% unname()
}


preprocess_simulation_results <- function(results) {
  results <- results %>%
    select(-c(submitted, started, done, error, mem.used, batch.id, log.file, job.hash, job.name, time.queued,
              time.running, problem, tags))
  
  betas = matrix(unlist(lapply(results$result, "[[", "beta")), ncol=4, byrow=TRUE) |> set_colnames(paste0("beta", 1:4))
  
  sigma = matrix(
    unlist(
      lapply(1:nrow(results), \(i) {
        if (results$algorithm[[i]] == "rankbased") {
          results$result[[i]][["sigma"]]
        } else {
          results$result[[i]][["variances"]][1]
        }
      })), ncol=1) |> set_colnames("sigma")
  
  theta = matrix(unlist(
    lapply(1:nrow(results), \(i) {
      if (results$algorithm[[i]] == "rankbased") {
        results$result[[i]][["theta"]]
      } else {
        results$result[[i]][["variances"]][2:3]
      }
    })), ncol=2, byrow=TRUE) |> set_colnames(paste0("theta", 0:1))
  
  init = matrix(
    unlist(
      lapply(1:nrow(results), function(i) { 
        if (results$algorithm[i] == "rankbased") { 
          results$result[[i]]$beta_init 
        } else { 
          rep(NA, 4)
        } 
      })), ncol=4, byrow=TRUE) |> set_colnames(paste0("beta", 1:4, "_init"))
  
  
  preprocessed  <- cbind(results, betas, sigma, theta, init)
  
  preprocessed <- preprocessed %>% select(-result) 
  
  preprocessed = preprocessed |>
    mutate(method = ifelse(algorithm != "rankbased", algorithm, ifelse(weighted, "weighted_rank", "rank"))) |>
    mutate(method = create_method_column(method))
  
  return(preprocessed)
}


colorcode <- c("REML" = "#458B74", "Rank" = "#FA681E", "SMDM" = "#0B4768",
               "Weighted Rank" = "#FF69B4")
ltycode <- c("REML" = "dotdash", "Rank" = "solid", "SMDM" = "dotted", "Weighted Rank" = "dashed")
shapecode <- c("REML" = 16, "Rank" = 25, "SMDM" = 3, "Weighted Rank" = 23) 

colors_efficiency_gains <- c("initial" = "gray", "updated" = "black")

my_palette <- RColorBrewer::brewer.pal(name="Greys",n=7)[2:7] %>%
  set_names(c(5, 10, 20, 30, 40, 50))

theme_set(theme_custom())


##### Functions for visualization of the results:

calc_mse <- function(x, true_value = 1) {
  
  mean((x - true_value)^2)
  
}


# Function to create lineplots as in the paper
lineplot_mean <- function(
    data, str,
    xaxis = "outlier_size",
    outlier_positioning = "random",
    includelme = TRUE,
    faceting_variable = "outlier_proportion",
    additional_filters = NULL, title = "", ylim = c(NA, NA),
    aggregation_function = "mean",
    variables = "coefficients",
    add_hline = TRUE,
    ylab = "Average estimate",
    xlab = "Outlier size",
    scales = "free"
) {
  
  if (variables == "coefficients") {
    
    dat <- data %>%
      pivot_longer(contains("beta") & !contains("init")) %>% 
      filter(if (isFALSE(includelme)) { method != "REML" } else { method != "bla" }) %>% 
      filter(if (!is.null(additional_filters)) { 
        rlang::eval_tidy(rlang::parse_expr(additional_filters)) 
      } else { method != "bla" }) %>%
      filter(single_outlier_type == outlier_positioning) %>%  
      # pivot_longer(all_of(paste0("beta", 1:4))) %>%
      group_by(outlier_size, outlier_level, outlier_proportion, method, outlier_type, name) %>%#
      summarize(value = do.call(aggregation_function, list(value))) %>% # mean((value - 1)^2)) %>% 
      mutate(coefficient = unlist(str_extract_all(name, "beta[1,2,3,4]"))) %>%
      mutate(coefficient = coefficients_tex(coefficient))
    
  } else if (variables == "variances") {
    
    dat <- data %>%
      pivot_longer(contains("sigma") | contains("theta")) %>% 
      mutate(true_value = ifelse(str_detect(name, "theta"), 0.5, 1)) %>%
      filter(if (!is.null(additional_filters)) { 
        rlang::eval_tidy(rlang::parse_expr(additional_filters)) 
      } else { method != "bla" })   %>%
      filter(if (isFALSE(includelme)) { method != "REML" } else { method != "bla" }) %>% 
      filter(single_outlier_type == outlier_positioning) %>%  
      # pivot_longer(all_of(paste0("beta", 1:4))) %>%
      group_by(outlier_size, outlier_level, outlier_proportion, method, outlier_type, name) %>%#
      summarize(value = do.call(aggregation_function, c(list(value), true_value = list(true_value)))) %>% # mean((value - 1)^2)) %>% 
      mutate(coefficient = unlist(str_extract_all(name, "theta[0, 1]|sigma"))) %>%
      mutate(coefficient = coefficients_tex(coefficient))
    
  }
  
  dat %>%
    #     rename("Outl. prop." = outlier_proportion) %>% 
    ggplot(aes_string(xaxis, "value", color = "method", linetype = "method", shape="method", fill="method")) +
    geom_line() +#alpha = 0.3) +
    geom_point() + #alpha = 0.3) +
    #  geom_smooth(se = FALSE) + 
    facet_wrap(. ~ coefficient, scales = scales, nrow = 1) +
    scale_color_manual(values = {if (includelme) colorcode else colorcode[-1]}) +
    scale_fill_manual(values = {if (includelme) colorcode else colorcode[-1]}) +
    scale_shape_manual(values = {if (includelme) shapecode else shapecode[-1]}) +
    scale_linetype_manual(values = {if (includelme) ltycode else ltycode[-1]}) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    ggtitle(title) +
    ylab(ylab) +
    xlab(xlab) +
    coord_cartesian(ylim = ylim) +
    { if (add_hline) {
      geom_hline(aes(yintercept = ifelse(str_detect(name, "theta"), 0.5, 1)), linetype = "dotted")
    } else NULL }
  
}


boxplots_consistency <- function(data, method, additional_filters=NULL) {
  (data %>%
     filter(method == method) |>
     filter(if (!is.null(additional_filters)) { 
       rlang::eval_tidy(rlang::parse_expr(additional_filters)) 
     } else { method != "bla" }) |>
     pivot_longer(contains("beta") & !contains("init")) %>%
     mutate(name = coefficients_tex(name)) %>%
     ggplot(aes(factor(n.groups), value, color = factor(group_size))) +
     geom_boxplot() +
     scale_color_manual(values = my_palette) +
     facet_wrap(name ~ ., nrow=1) +#, scales = "free") +
     geom_hline(aes(yintercept = ifelse(str_detect(name, "theta"), 0.5, 1)), linetype = "dotted") +
     xlab("g") + ylab("Coefficient estimates") +
     guides(color = guide_legend(nrow = 1, title = "$n_i$:")) +
     ylim(c(0, 2.5))
   ) +
    (data %>%
       filter(method == method) |>
       filter(if (!is.null(additional_filters)) { 
         rlang::eval_tidy(rlang::parse_expr(additional_filters)) 
       } else { method != "bla" }) |> 
       pivot_longer(c(contains("sigma"), contains("theta"))) %>%
       mutate(name = coefficients_tex(name)) %>%
       mutate(name = factor(name, 
                            levels = c("$\\theta_0$", "$\\theta_1$", "$\\sigma$"))) %>% 
       ggplot(aes(factor(n.groups), value, color = factor(group_size))) +
       geom_boxplot() +
       facet_wrap(name ~ ., scales = "free", nrow = 1) +
       geom_hline(aes(yintercept = ifelse(str_detect(name, "theta"), 0.5, 1)), linetype = "dotted") +
       scale_color_manual(values = my_palette) +
       xlab("g") + ylab("Coefficient estimates") +
       guides(color = guide_legend(nrow = 1, title = "$n_i$:")))  +
    plot_layout(guides = "collect", nrow = 2)
}


plot_mse_efficiency <- function(dat, title, ylim = c(NA, NA)) {
  dat %>%
    filter(method == "Rank") %>%
    pivot_longer(contains("beta")) |>
    mutate(coefficient = unlist(coefficients_tex(name)),
           initial = ifelse(str_detect(name, "init"), "initial", "updated")) %>%
    select(-name) %>%
    group_by(coefficient, initial, outlier_size, outlier_proportion) %>%
    summarize(value = mean((value - 1)^2)) %>%
    ggplot(aes(outlier_size, value, color = initial)) +
    geom_line(alpha = 0.5) +
    geom_point(alpha = 0.5) + 
    geom_smooth(se = FALSE) +
    facet_wrap(coefficient ~ ., nrow = 1) +
    scale_color_manual(values = colors_efficiency_gains) +
    theme(legend.title = element_blank()) +
    xlab("Outlier size") +
    ylab("MSE") +
    ggtitle(title) +
    coord_cartesian(ylim = ylim)
}


