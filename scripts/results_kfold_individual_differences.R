rm(list = ls())

# Notes -------------------------------------------------------------------

  #########################################################################
  # Adjust n_cores to the number of CPU cores you want to use.            #
  # WARNING: Running this script might take several hours/days depending  #
  #          on the computer you use.                                     #
  #########################################################################

# Library -----------------------------------------------------------------

  # Load packages
  library(tidyverse); library(rstan); library(loo)
  
  # Stan options
  n_cores <- 8L
  options(mc.cores = n_cores)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
  rstan_options(auto_write = FALSE)
  seeds <- c(5477649, 3984665, 8989566, 9363790)
  

# Functions ---------------------------------------------------------------
  inv_logit <- function(x) exp(x) / ( 1 + exp(x) )
  
  
# Prepare -----------------------------------------------------------------

  # Import data
  dl <- read_rds("data/d1.rds") %>%
    filter(!is.na(q1)) %>%
    mutate(og = case_when(
      category == 1 & ig_caste == "gm" ~ 0L,
      category == 2 & ig_caste == "obc" ~ 0L,
      category == 3 & ig_caste == "scst" ~ 0L,
      TRUE ~ 1L
    )) %>%
    arrange(participant, target)

  # Import imputed data
  dl_imp <- read_rds("results/ic_impute.rds") %>%
    mutate(og = case_when(
      category == 1 & ig_caste == "gm" ~ 0L,
      category == 2 & ig_caste == "obc" ~ 0L,
      category == 3 & ig_caste == "scst" ~ 0L,
      TRUE ~ 1L
    )) %>%
    select(-z_obs, -z_imp)

  # Summarize/Standardize
  dl_imp <- dl_imp %>%
    filter(!is.na(x_obs)) %>%
    select(-.draw) %>%
    crossing(.draw = 1:max(dl_imp$.draw, na.rm = TRUE)) %>%
    bind_rows(dl_imp %>% filter(is.na(x_obs))) %>%
    group_by_at(vars(-x_obs, -x_imp)) %>%
    summarise_at(vars(x_obs, x_imp), mean) %>%
    group_by(.draw, item) %>%
    mutate(
      x_obs = if_else(og == 0L, NA_real_, x_obs),
      z_obs = (x_obs - mean(x_obs, na.rm = TRUE))/sd(x_obs, na.rm = TRUE),
      z_imp = (x_imp - mean(x_obs, na.rm = TRUE))/sd(x_obs, na.rm = TRUE)
    ) %>%
    ungroup()

  # Summarize posterior distributions for model
  dl_imp <- dl_imp %>%
    group_by(participant, category, item) %>%
    summarize(
      x_obs = unique(x_obs),
      z_obs = unique(z_obs),
      z_mu = mean(z_imp),
      z_sigma = sd(z_imp)
    ) %>%
    ungroup()

  # Merge missing and observed data
  dl_imp <- dl_imp %>%
    complete(participant, item, category = 1:6) %>%
    left_join(
      dl %>% select(participant, target, category, og),
      .,
      by = c("participant", "category")
    ) %>%
    mutate(
      z_obs = if_else(og == 0L | category > 4L, -99, z_obs)
    )

  # Arrange for model
  dl_imp <- dl_imp %>%
    arrange(item, participant, target) %>%
    rowid_to_column("ii")

  # Add factor scores for SDO-D and SDO-E
  source("scripts/prepare_data_sdo.R")

  # Compose data list
  data_list <- with(dl, list(
      N = length(participant),
      J = n_distinct(participant),
      K = n_distinct(category),
      jj = participant,
      kk = category,
      x_scst = ifelse(ig_caste == "scst", 1L, 0L),
      x_obc = ifelse(ig_caste == "obc", 1L, 0L),
      x_og = og,
      x_gs = gs,
      x_sdo_d = sdo_d,
      x_sdo_e = sdo_e,
      y = q1
    )) %>%
    append(
      with(dl_imp, list(
        N_obs = sum(!is.na(z_obs)),
        N_mis = sum(is.na(z_obs)),
        ii_obs = ii[!is.na(z_obs)],
        ii_mis = ii[is.na(z_obs)],
        x_obs = z_obs[!is.na(z_obs)],
        x_mu = z_mu[is.na(z_obs)],
        x_sigma = z_sigma[is.na(z_obs)]
      ))
    )
  
  # Export/Import
  # write_rds(data_list, "dlist_for_kfold.rds")
  # data_list <- read_rds("dlist_for_kfold.rds")

  # Generate indeces for 10-fold cross-validation
  set.seed(24987504)
  d_kfold <- tibble(
    ii   = 1:data_list$N, 
    jj   = data_list$jj, 
    fold = kfold_split_stratified(K = 10, x = jj)
  )

  
# Functions ---------------------------------------------------------------
  
  # Append indeces for 10-fold cross-validation
  data_fold <- function(k) {
    append(data_list, list(
      N_t = d_kfold %>% filter(fold != k) %>% nrow(),    # training
      N_h = d_kfold %>% filter(fold == k) %>% nrow(),    # hold out
      ii_t = d_kfold %>% filter(fold != k) %>% pull(ii), # training
      ii_h = d_kfold %>% filter(fold == k) %>% pull(ii)  # hold out
    ))
  }


# Run ---------------------------------------------------------------------

  # Model: M2 + Participant: Intergroup contact (4)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m4_kfold.stan")
    print(paste0("M4: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 1000 + 8000/n_cores, warmup = 1000,
      chains = n_cores,
      seed = seeds[1],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m4_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m4_llk <- append(m4_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m4_llk, "models/kfold/q1_m4_llk.rds")


  # Model: M2 + Participant: Intergroup contact (2)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m5_kfold.stan")
    print(paste0("M5: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 1000 + 8000/n_cores, warmup = 1000,
      chains = n_cores,
      seed = seeds[2],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m5_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m5_llk <- append(m5_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m5_llk, "models/kfold/q1_m5_llk.rds")

  # Model: M5 + Participant: Intergroup contact (2) * Target: Category (4)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m6_kfold.stan")
    print(paste0("M6: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 1000 + 8000/n_cores, warmup = 1000,
      chains = n_cores,
      control = list(adapt_delta = 0.95),
      seed = seeds[3],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m6_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m6_llk <- append(m6_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m6_llk, "models/kfold/q1_m6_llk.rds")
  
  # Model: M2 + Participant: SDO (2)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m7_kfold.stan")
    print(paste0("M7: ", k, "/10"))
    fit <- sampling(
      model, 
      data = data_fold(k), 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      seed = seeds[4],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m7_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m7_llk <- append(m7_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m7_llk, "models/kfold/q1_m7_llk.rds")
  

# Extract -----------------------------------------------------------------

  # Load model results
  m4_llk <- read_rds("models/kfold/q1_m4_llk.rds")
  m5_llk <- read_rds("models/kfold/q1_m5_llk.rds")
  m6_llk <- read_rds("models/kfold/q1_m6_llk.rds")
  m7_llk <- read_rds("models/kfold/q1_m7_llk.rds")
  
  # Compile model results
  q1_kfold <- list(m4_llk, m5_llk, m6_llk, m7_llk)

  # Merge with m0_llk to m3_llk
  q1_kfold <- append(read_rds("results/q1_kfold.rds"), q1_kfold)
  
  # Export
  write_rds(q1_kfold, "results/q1_kfold.rds")

  # Calculate expected log pointwise predictive density for each model
  q1_elpd <- q1_kfold %>% 
    map(., ~map(., ~log(colMeans(exp(.))))) %>%
    map(., unlist) %>% 
    tibble(
      model = 0:(length(q1_kfold)-1),
      elpd_i = .
    )
  
  # Export
  write_rds(q1_elpd, "results/q1_elpd.rds")
  