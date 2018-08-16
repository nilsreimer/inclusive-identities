rm(list = ls())

# Notes -------------------------------------------------------------------
  # Explain that you use single imputation (best guess) for intergroup contact
  # because missingness is negligible (n = 66 out of ...) and explain how it doesn't make a difference.

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes); library(rstan); library(loo)

  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seed = 0892561437


# Import ------------------------------------------------------------------
  source("prepare_data_intergroup_contact.R")

  # Compose data list
  data_list <- list(
    # Numbers
    N_row = nrow(dl),
    N_col = n_distinct(dl_imp$item),
    N_mis = sum(dl_imp$missing),
    N_obs = sum(!dl_imp$missing),
    J = n_distinct(dl$participant),
    K = n_distinct(dl$category),
    # Indices
    ii_mis  = dl_imp$ii[dl_imp$missing],
    ii_obs  = dl_imp$ii[!dl_imp$missing],
    jj = dl$participant,
    kk = dl$category,
    # Vectors
    x_scst = ifelse(dl$ig_caste == "scst", 1L, 0L),
    x_obc = ifelse(dl$ig_caste == "obc", 1L, 0L),
    x_og = dl$og,
    x_imp_m  = dl_imp$x_imp[dl_imp$missing],
    x_imp_sd = dl_imp$x_imp_sd[dl_imp$missing], 
    x_obs = dl_imp$x_imp[!dl_imp$missing],
    y = dl$q1
  )


# Models ------------------------------------------------------------------
  # Model: 
  m4_stan <- stan("models/q1_m4.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000,
                  pars = c("alpha", "x", "X", "x_imp"), include = FALSE)
  
  # Model: 
  m5_stan <- stan("models/q1_m5.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000,
                  pars = c("alpha", "x", "X", "x_imp"), include = FALSE)
  
  # Model: 
  m6_stan <- stan("models/q1_m6.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000,
                  pars = c("alpha", "x", "X", "x_imp"), include = FALSE)
  
  m4_stan %>% 
    gather_draws(b_cq, b_pc, b_nc, b_of) %>%
    mutate(.value = exp(.value)) %>%
    median_hdi(.width = .97)

  write_rds(m4_stan, "models/m4_stan_imputed.rds")

# Compare -----------------------------------------------------------------
  # Compare models
  mc_fit <- list(m0_stan, m1_stan, m2_stan, m3_stan)
  mc_llk <- map(mc_fit, extract_log_lik, merge_chains = FALSE)
  mc_loo <- map(mc_llk, ~loo(., r_eff = relative_eff(exp(.)), cores = 4))

  compare(x = mc_loo)
  
  # Export
  write_rds(mc_fit, "models/q1_mc_fit.rds")
  write_rds(mc_loo, "models/q1_mc_loo.rds")


# Predict -----------------------------------------------------------------
  # Extract posterior samples
  post <- tibble(
    model = 0:3,
    fit = mc_fit
  ) %>% mutate(
    p_est = map(fit, ~spread_draws(., p_k[category, ig_caste]))
  ) %>% select(model, p_est)
  
  # Export
  write_rds(post, "models/q1_post.rds")
