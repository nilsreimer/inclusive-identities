rm(list = ls())

# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(loo)
  
  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seeds = c(5477649, 3984665, 8989566, 9363790)


# Import ------------------------------------------------------------------
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
    x_scst = ifelse(ig_caste == "scst", 1L, 0L),
    x_obc = ifelse(ig_caste == "obc", 1L, 0L),
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
  # Model 0: Varying Intercept (Participant)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m0_kfold.stan",
      data = data_fold(k), 
      seed = seeds[1],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m0_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m0_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m0_llk_t <- append(m0_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m0_llk_h <- append(m0_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m0_llk_h, "models/kfold/q1_m0_llk_h.rds")
  
  # Model 1: + Varying Intercept (Target: Category)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m1_kfold.stan",
      data = data_fold(k),
      seed = seeds[2],
      iter = 1000,
      warmup = 750,
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m1_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m1_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m1_llk_t <- append(m1_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m1_llk_h <- append(m1_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m1_llk_h, "models/kfold/q1_m1_llk_h.rds")

  # Model 2: + Participant: SCST * Target: Category (4)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m2_kfold.stan",
      data = data_fold(k),
      seed = seeds[3],
      iter = 1000,
      warmup = 750,
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m2_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m2_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m2_llk_t <- append(m2_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m2_llk_h <- append(m2_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m2_llk_h, "models/kfold/q1_m2_llk_h.rds")
  
  # Model 3: + Participant: OBC * Target: Category (4)
  for (k in 1:10) {
    print(paste("[", k, "/10]", sep = ""))
    fit <- stan(
      "models/kfold/q1_m3_kfold.stan",
      data = data_fold(k), 
      seed = seeds[4],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h", "log_lik_t")
    )
    if (k == 1) {
      m3_llk_t <- list(extract_log_lik(fit, parameter_name = "log_lik_t"))
      m3_llk_h <- list(extract_log_lik(fit, parameter_name = "log_lik_h"))
    } else {
      m3_llk_t <- append(m3_llk_t, list(extract_log_lik(fit, parameter_name = "log_lik_t")))
      m3_llk_h <- append(m3_llk_h, list(extract_log_lik(fit, parameter_name = "log_lik_h")))
    }
  }
  write_rds(m3_llk_h, "models/kfold/q1_m3_llk_h.rds") 


# Estimate ----------------------------------------------------------------
  m0_llk_h <- read_rds("models/kfold/q1_m0_llk_h.rds")
  m1_llk_h <- read_rds("models/kfold/q1_m1_llk_h.rds")
  m2_llk_h <- read_rds("models/kfold/q1_m2_llk_h.rds")
  m3_llk_h <- read_rds("models/kfold/q1_m3_llk_h.rds")
  
  # Calculate expected log pointwise predictive density for each model
  llk_h  <- list(m0_llk_h, m1_llk_h, m2_llk_h, m3_llk_h)
  elpd_i <- map(llk_h, ~map(., ~log(colMeans(exp(.))))) %>%
            map(., unlist)
  
  # Compare models
  elpd_kfold <- tibble(
    model  = 1:length(elpd_i),
    elpd_i = elpd_i
  ) %>% 
  mutate(
    elpd    = map(elpd_i, sum) %>% unlist(),
    elpd_se = map(elpd_i, ~(sd(.) * sqrt(length(.)))) %>% unlist(),
    elpd_d  = map2(lag(elpd_i), elpd_i, compare_kfold)
  ) %>% unnest(elpd_d)

  # Export
  elpd_kfold %>% write_rds(., "models/kfold/q1_elpd_kfold.rds")
