rm(list = ls())

# Notes -------------------------------------------------------------------

  #########################################################################
  # Adjust n_cores to the number of CPU cores you want to use.            #
  # For this analysis, we are only using observed responses.              #
  # WARNING: Running this script might take several hours/days depending  #
  #          on the computer you use.                                     #
  #########################################################################

# Library -----------------------------------------------------------------
  
  # Load packages
  library(tidyverse); library(rstan); library(tidybayes); library(loo)

  # Stan options
  n_cores <- 8L
  options(mc.cores = n_cores)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
  rstan_options(auto_write = TRUE)
  seeds = c(1157744, 2691555, 4045395, 2910499, 5195618, 5975015, 2391421)


# Import ------------------------------------------------------------------
  
  # Import data
  dl <- read_rds("data/d1.rds") %>%
    filter(!is.na(q1)) %>%
    mutate_at(
      vars(q2, q3), 
      list("mean" = mean, "sd" = sd), 
      na.rm = TRUE
    ) %>%
    mutate(
      q2_z = (q2 - q2_mean) / q2_sd, 
      q3_z = (q3 - q3_mean) / q3_sd
    ) %>%
    arrange(participant, target)
  
  # Divides data for 10-fold cross-validation
  set.seed(24987504)
  dl <- dl %>% mutate(fold = kfold_split_stratified(K = 10, x = participant))

  # Exclude cases with missing responses on either outcome variable
  dl <- dl %>% filter(!is.na(q2), !is.na(q3))

  # Compose data list
  data_list <- with(dl, list(
      N_row = length(participant),
      N_col = 2L,
      J = n_distinct(participant),
      K = n_distinct(category),
      jj = as.integer(factor(participant)),
      kk = category,
      x_scst = if_else(ig_caste == "scst", 1L, 0L),
      x_obc = if_else(ig_caste == "obc", 1L, 0L),
      x_q1 = q1,
      Y = cbind(q2_z, q3_z)
    ))
  
  
# Functions ---------------------------------------------------------------
  
  # Append indeces for 10-fold cross-validation
  data_fold <- function(k) {
    append(data_list, with(
      dl %>% transmute(ii = row_number(), fold), 
      list(
        N_t = length(ii[fold != k]), # training
        N_h = length(ii[fold == k]), # hold out
        ii_t = ii[fold != k],        # training
        ii_h = ii[fold == k]         # hold out
      )
    ))
  }


# Models ------------------------------------------------------------------
  
  # Model: Varying intercept (Participant)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q2q3_m0_kfold.stan")
    print(paste0("M0: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 750 + 2000/n_cores, warmup = 750,
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
  write_rds(m0_llk, "models/kfold/q2q3_m0_llk.rds")
  m0_llk <- read_rds("models/kfold/q2q3_m0_llk.rds")

  # Model: M0 + Varying Intercept (Target: Category)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q2q3_m1_kfold.stan")
    print(paste0("M1: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 750 + 2000/n_cores, warmup = 750,
      control = list(max_treedepth = 12),
      chains = n_cores,
      seed = seeds[2],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m1_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m1_llk <- append(m1_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m1_llk, "models/kfold/q2q3_m1_llk.rds")

  # Model: M1 + Participant: SCST * Target: Category (4)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q2q3_m2_kfold.stan")
    print(paste0("M2: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 750 + 2000/n_cores, warmup = 750,
      chains = n_cores,
      control = list(max_treedepth = 12),
      seed = seeds[3],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m2_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m2_llk <- append(m2_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m2_llk, "models/kfold/q2q3_m2_llk.rds")

  # Model: M2 + Participant: OBC * Target: Category (4)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q2q3_m3_kfold.stan")
    print(paste0("M3: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 750 + 2000/n_cores, warmup = 750,
      chains = n_cores,
      control = list(max_treedepth = 12),
      seed = seeds[4],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m3_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m3_llk <- append(m3_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m3_llk, "models/kfold/q2q3_m3_llk.rds")

  # Model: M3 + Fixed effect (Categorization)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q2q3_m4_kfold.stan")
    print(paste0("M4: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 750 + 2000/n_cores, warmup = 750,
      chains = n_cores,
      control = list(max_treedepth = 12),
      seed = seeds[5],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m4_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m4_llk <- append(m4_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m4_llk, "models/kfold/q2q3_m4_llk.rds")

  # Model: M4 + Varying effect (Categories)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q2q3_m5_kfold.stan")
    print(paste0("M5: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 750 + 2000/n_cores, warmup = 750,
      chains = n_cores,
      control = list(max_treedepth = 12),
      seed = seeds[6],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m5_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m5_llk <- append(m5_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m5_llk, "models/kfold/q2q3_m5_llk.rds")

  # Model: M5 + Interaction (SC/ST, OBC)
  for (k in 1:10) {
    if (k == 1) model <- stan_model("models/kfold/q2q3_m6_kfold.stan")
    print(paste0("M6: ", k, "/10"))
    fit <- sampling(
      model,
      data = data_fold(k),
      iter = 750 + 2000/n_cores, warmup = 750,
      chains = n_cores,
      control = list(max_treedepth = 12),
      seed = seeds[7],
      pars = c("log_lik"), include = TRUE
    )
    if (k == 1) {
      m6_llk <- list(extract_log_lik(fit, parameter_name = "log_lik"))
    } else {
      m6_llk <- append(m6_llk, list(extract_log_lik(fit, parameter_name = "log_lik")))
    }
  }
  write_rds(m6_llk, "models/kfold/q2q3_m6_llk.rds")
  
# Extract -----------------------------------------------------------------

  # Load model results
  m0_llk <- read_rds("models/kfold/q2q3_m0_llk.rds")
  m1_llk <- read_rds("models/kfold/q2q3_m1_llk.rds")
  m2_llk <- read_rds("models/kfold/q2q3_m2_llk.rds")
  m3_llk <- read_rds("models/kfold/q2q3_m3_llk.rds")
  m4_llk <- read_rds("models/kfold/q2q3_m4_llk.rds")
  m5_llk <- read_rds("models/kfold/q2q3_m5_llk.rds")
  m6_llk <- read_rds("models/kfold/q2q3_m6_llk.rds")
  
  # Compile model results
  q2q3_kfold <- list(m0_llk, m1_llk, m2_llk, m3_llk, m4_llk, m5_llk, m6_llk)

  # Export
  write_rds(q2q3_kfold, "results/q2q3_kfold.rds")

  # Calculate expected log pointwise predictive density for each model
  q2q3_elpd <- q2q3_kfold %>%
    map(., ~map(., ~log(colMeans(exp(.))))) %>%
    map(., unlist) %>%
    tibble(
      model = 0:(length(q2q3_kfold)-1),
      elpd_i = .
    )

  # Export
  write_rds(q2q3_elpd, "results/q2q3_elpd.rds")
