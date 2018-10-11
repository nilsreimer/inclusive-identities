rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(loo)
  
  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seeds = c(6678002, 8970252,	9556204, 8660009)


# Import ------------------------------------------------------------------
  source("data/prepare_data_ic.R")
  source("data/prepare_data_sdo.R")

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
    x_gs = gs,
    x_sdo_d = sdo_d,
    x_sdo_e = sdo_e,
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
  data_fold <- function(k) {
    append(data_list, list(
      N_t = d_kfold %>% filter(fold != k) %>% nrow(),
      N_h = d_kfold %>% filter(fold == k) %>% nrow(),
      ii_t = d_kfold %>% filter(fold != k) %>% pull(ii),
      ii_h = d_kfold %>% filter(fold == k) %>% pull(ii)
    ))
  }
  
  compare_kfold <- function(elpd1, elpd2) {
    d_elpd  <- sum(elpd2 - elpd1)
    d_se    <- sd(elpd2 - elpd1) * sqrt(length(elpd1))
    d_ratio <- d_elpd / d_se
    return(tibble(d_elpd = d_elpd, d_se = d_se, d_ratio = d_ratio))
  }
  
  
# Models ------------------------------------------------------------------
  # Model 4: M2 + Participant: Intergroup contact (4)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m4_kfold.stan",
      data = data_fold(k), 
      seed = seeds[1],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m4_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m4_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m4_llk_t <- append(m4_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m4_llk_h <- append(m4_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m4_llk_h, "models/kfold/q1_m4_llk_h.rds") 
  
  # Model 5: M2 + Participant: Intergroup contact (2)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m5_kfold.stan",
      data = data_fold(k), 
      seed = seeds[2],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m5_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m5_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m5_llk_t <- append(m5_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m5_llk_h <- append(m5_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m5_llk_h, "models/kfold/q1_m5_llk_h.rds") 
  
  # Model: M6 + Participant: Intergroup contact (2) * Target: Category (4)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m6_kfold.stan",
      data = data_fold(k), 
      seed = seeds[3],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m6_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m6_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m6_llk_t <- append(m6_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m6_llk_h <- append(m6_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m6_llk_h, "models/kfold/q1_m6_llk_h.rds") 
  
  # Model 7: M2 + Participant: SDO (2)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m7_kfold.stan",
      data = data_fold(k), 
      seed = seeds[4],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m7_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m7_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m7_llk_t <- append(m7_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m7_llk_h <- append(m7_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m7_llk_h, "models/kfold/q1_m7_llk_h.rds") 


# Estimate ----------------------------------------------------------------
  m0_llk_h <- read_rds("models/kfold/q1_m0_llk_h.rds")
  m1_llk_h <- read_rds("models/kfold/q1_m1_llk_h.rds")
  m2_llk_h <- read_rds("models/kfold/q1_m2_llk_h.rds")
  m3_llk_h <- read_rds("models/kfold/q1_m3_llk_h.rds")
  m4_llk_h <- read_rds("models/kfold/q1_m4_llk_h.rds")
  m5_llk_h <- read_rds("models/kfold/q1_m5_llk_h.rds")
  m6_llk_h <- read_rds("models/kfold/q1_m6_llk_h.rds")
  m7_llk_h <- read_rds("models/kfold/q1_m7_llk_h.rds")
  
  # Calculate expected log pointwise predictive density for each model
  llk_h  <- list(m0_llk_h, m1_llk_h, m2_llk_h, m3_llk_h, 
                 m4_llk_h, m5_llk_h, m6_llk_h, m7_llk_h)
  elpd_i <- map(llk_h, ~map(., ~log(colMeans(exp(.))))) %>%
            map(., unlist)
  
  # Compare models
  elpd_kfold <- tibble(
    model  = 0:(length(elpd_i)-1),
    elpd_i = elpd_i
  ) %>% 
    mutate(
      elpd    = map(elpd_i, sum) %>% unlist(),
      elpd_se = map(elpd_i, ~(sd(.) * sqrt(length(.)))) %>% unlist(),
      elpd_d  = case_when(
        model == 4 ~ map2(lag(elpd_i, 2), elpd_i, compare_kfold),
        model == 7 ~ map2(lag(elpd_i, 5), elpd_i, compare_kfold),
        TRUE ~ map2(lag(elpd_i), elpd_i, compare_kfold)
      )
    ) %>% unnest(elpd_d)
  
  # Export
  elpd_kfold %>% write_rds(., "models/kfold/q1_elpd_kfold.rds")
  
  # Calculate stacking weights
  matrix(unlist(elpd_i), ncol = length(elpd_i)) %>% stacking_weights()
  matrix(unlist(elpd_i), ncol = length(elpd_i)) %>% pseudobma_weights()
