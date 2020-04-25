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


# Prepare -----------------------------------------------------------------

  # Import data
  dl <- read_rds("data/d1.rds") %>%
        filter(!is.na(q1)) %>%
        arrange(participant, target)
  
  # Compose data list
  data_list <- with(dl, list(
    N = length(participant),
    J = n_distinct(participant),
    K = n_distinct(category),
    jj = participant,
    kk = category,
    x_scst = if_else(ig_caste == "scst", 1L, 0L),
    x_obc = if_else(ig_caste == "obc", 1L, 0L),
    y = q1
  ))
  
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

  # Model: Varying Intercept (Participant)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m0_kfold.stan")
    print(paste0("M0: ", k, "/10"))
    fit <- sampling(
      model, 
      data = data_fold(k), 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      seed = seeds[1],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m0_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m0_llk <- append(m0_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m0_llk, "models/kfold/q1_m0_llk.rds")
  
  # Model: + Varying Intercept (Target: Category)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m1_kfold.stan")
    print(paste0("M1: ", k, "/10"))
    fit <- sampling(
      model, 
      data = data_fold(k), 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      control = list(adapt_delta = 0.90),
      seed = seeds[2],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m1_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m1_llk <- append(m1_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m1_llk, "models/kfold/q1_m1_llk.rds")
  
  # Model: + Participant: SCST * Target: Category (4)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m2_kfold.stan")
    print(paste0("M2: ", k, "/10"))
    fit <- sampling(
      model, 
      data = data_fold(k), 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      seed = seeds[3],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m2_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m2_llk <- append(m2_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m2_llk, "models/kfold/q1_m2_llk.rds")
  
  # Model: + Participant: OBC * Target: Category (4)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q1_m3_kfold.stan")
    print(paste0("M3: ", k, "/10"))
    fit <- sampling(
      model, 
      data = data_fold(k), 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      control = list(adapt_delta = 0.95),
      seed = seeds[4],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m3_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m3_llk <- append(m3_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m3_llk, "models/kfold/q1_m3_llk.rds")


# Extract -----------------------------------------------------------------

  # Compile model results
  q1_kfold <- list(m0_llk, m1_llk, m2_llk, m3_llk)
  
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
  write_rds(q1_elpd, "results/q1_elpd")
