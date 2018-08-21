rm(list = ls())

# Notes -------------------------------------------------------------------
  # In total, I considered 5 * 4 = 20 items measuring intergroup contact. Across
  # items, less than 1.6% (n = 66) of individual responses were missing - though
  # about 20% of all participants failed to respond to at least one item. As 
  # item-by-item missingness was negligible, I used single imputation rather 
  # than (computationally demanding) full imputation. Before running Models 4 to 
  # 7, I imputed missing values from a multivariate normal distribution (see 
  # "impute_data.R"). For all subsequent analyses, I substituted the most likely
  # value (from the posterior distribution of the imputation model) for each 
  # missing value (in the original dataset).

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes); library(rstan); library(loo)

  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seed = 0892561437


# Import ------------------------------------------------------------------
  source("data/prepare_data_ic.R")

  # Compose data list
  data_list <- with(dl, list(
    # Numbers
    N = length(participant),
    J = n_distinct(participant),
    K = n_distinct(category),
    # Indices
    jj = participant,
    kk = category,
    # Vectors
    x_scst = ifelse(ig_caste == "scst", 1L, 0L),
    x_obc = ifelse(ig_caste == "obc", 1L, 0L),
    x_og = og,
    x_cq = replace_na(cq, -99), # Replace NAs for foreign targets 
    x_pc = replace_na(pc, -99), # Replace NAs for foreign targets
    x_nc = replace_na(nc, -99), # Replace NAs for foreign targets
    x_of = replace_na(of, -99), # Replace NAs for foreign targets
    y = q1
  ))


# Models ------------------------------------------------------------------
  # Model: M2 + Participant: Intergroup contact (4)
  m4_stan <- stan("models/q1_m4.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Model: M2 + Participant: Intergroup contact (2)
  m5_stan <- stan("models/q1_m5.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Model: M5 + Participant: Intergroup contact (2) * Target: Category (4)
  m6_stan <- stan("models/q1_m6.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  m4_stan %>% 
    gather_draws(b_cq, b_pc, b_nc, b_of) %>%
    mutate(.value = exp(.value)) %>%
    median_hdi(.width = .97)


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
