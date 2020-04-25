rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------

  # Load packages
  library(tidyverse); library(loo)


# Prepare -----------------------------------------------------------------

  # Import results from k-fold cross-validation
  q2q3_elpd <- read_rds("results/q2q3_elpd.rds") %>% filter(model %in% 0:6)
  
  # Calculate pseudo-BMA weights
  q2q3_elpd <- q2q3_elpd %>% 
    unnest(elpd_i) %>% 
    group_by(model) %>% 
    mutate(ii = row_number()) %>% 
    ungroup() %>% 
    pivot_wider(
      names_from = model, 
      names_prefix = "M", 
      values_from = elpd_i
    ) %>% 
    select(-ii) %>% 
    as.matrix() %>% 
    pseudobma_weights() %>% 
    mutate(q2q3_elpd, elpd_w = .)
  
  # Calculate expected log predictive density (ELPD)
  q2q3_elpd <- q2q3_elpd %>% 
    mutate(
      elpd    = map(elpd_i, sum) %>% unlist(),
      elpd_se = map(elpd_i, ~(sd(.) * sqrt(length(.)))) %>% unlist()
    )
  
  # Calculate differences in difference in ELPD
  q2q3_elpd <- q2q3_elpd %>% 
    select(model, elpd_i) %>% 
    crossing(comparison_model = 0:6) %>% 
    left_join(
      q2q3_elpd %>% select(comparison_model = model, comparison_elpd_i = elpd_i),
      by = "comparison_model"
    ) %>% 
    mutate(
      elpd_d = map2_dbl(elpd_i, comparison_elpd_i, ~sum(.x - .y)),
      elpd_d_se = map2_dbl(elpd_i, comparison_elpd_i, ~(sd(.x - .y) * sqrt(length(.x)))),
      elpd_d_z = if_else(elpd_d == 0, 0, elpd_d/elpd_d_se)
    ) %>% 
    select(model, comparison_model, elpd_d_z) %>% 
    pivot_wider(
      names_from = comparison_model, 
      names_prefix = "M", 
      values_from = elpd_d_z
    ) %>% 
    left_join(q2q3_elpd, ., by = "model")
  
  # Reformat
  q2q3_elpd <- q2q3_elpd %>% mutate(model = paste0("M", model))
  
  # View
  q2q3_elpd %>% 
    select(model, elpd_w, matches("M[0-9]")) %>% 
    mutate_at(vars(matches("M[0-9]")), round, digits = 1) %>% 
    mutate_at(vars(elpd_w), round, digits = 2)
