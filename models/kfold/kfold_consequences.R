rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(tidybayes); library(loo)
  
  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seeds = c(1157744, 2691555, 4045395, 2910499, 5195618, 5975015, 2391421)


# Import ------------------------------------------------------------------
  dl <- read_rds("data/d1.rds") %>%
        filter(!is.na(q1)) %>%
        mutate_at(vars(q2, q3), funs(mean, sd), na.rm = TRUE) %>%
        mutate(q2_z = (q2 - q2_mean) / q2_sd, q3_z = (q3 - q3_mean) / q3_sd) %>%
        arrange(participant, target)
  
  # Single (median) imputation
  Y           <- as.matrix(dl[,c("q2_z", "q3_z")])
  Y[is.na(Y)] <- read_rds("models/q2q3_mc_fit.rds")[[6]] %>%
                 spread_draws(y_mis[ii]) %>%
                 summarise(y_mis = median(y_mis)) %>%
                 pull(y_mis)

  # Compose data list
  data_list <- with(dl, list(
    # Numbers
    N_row = length(participant),
    N_col = ncol(Y),
    J = n_distinct(participant),
    K = n_distinct(category),
    # Indices
    jj = participant,
    kk = category,
    # Vectors
    x_scst = ifelse(ig_caste == "scst", 1L, 0L),
    x_obc = ifelse(ig_caste == "obc", 1L, 0L),
    x_q1 = q1,
    Y = Y
  ))
  
  # Generate indeces for 10-fold cross-validation
  set.seed(24987504)
  d_kfold <- tibble(
    ii   = 1:data_list$N_row, 
    jj   = data_list$jj, 
    fold = kfold_split_stratified(K = 10, x = jj),
    na   = is.na(dl$q2_z) | is.na(dl$q3_z)
  ) %>% mutate(
    fold = ifelse(na, 0, fold)
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
  # # Model 0
  # for (k in 1:10) {
  #   print(paste("[", k, "/10] - Model 1", sep = ""))
  #   fit <- stan(
  #     "models/kfold/q2q3_m0_kfold.stan",
  #     data = data_fold(k),
  #     seed = seeds[1],
  #     iter = 1000,
  #     warmup = 750,
  #     pars = c("log_lik_h")
  #   )
  #   if (k == 1) {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>%
  #     list(.) %>%
  #     write_rds(., "models/kfold/q2q3_m0_llk.rds")
  #   } else {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>%
  #     list(.) %>%
  #     append(read_rds("models/kfold/q2q3_m0_llk.rds"), .) %>%
  #     write_rds(., "models/kfold/q2q3_m0_llk.rds")
  #   }
  # }

  # # Model 1
  # for (k in 1:10) {
  #   print(paste("[", k, "/10] - Model 1", sep = ""))
  #   fit <- stan(
  #     "models/kfold/q2q3_m1_kfold.stan",
  #     data = data_fold(k),
  #     seed = seeds[2],
  #     iter = 1000,
  #     warmup = 750,
  #     pars = c("log_lik_h")
  #   )
  #   if (k == 1) {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>%
  #     list(.) %>%
  #     write_rds(., "models/kfold/q2q3_m1_llk.rds")
  #   } else {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>%
  #     list(.) %>%
  #     append(read_rds("models/kfold/q2q3_m1_llk.rds"), .) %>%
  #     write_rds(., "models/kfold/q2q3_m1_llk.rds")
  #   }
  # }

  # # Model 2
  # for (k in 1:10) {
  #   print(paste("[", k, "/10] - Model 2", sep = ""))
  #   fit <- stan(
  #     "models/kfold/q2q3_m2_kfold.stan",
  #     data = data_fold(k),
  #     seed = seeds[3],
  #     iter = 1000,
  #     warmup = 750,
  #     pars = c("log_lik_h")
  #   )
  #   if (k == 1) {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>%
  #     list(.) %>%
  #     write_rds(., "models/kfold/q2q3_m2_llk.rds")
  #   } else {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>%
  #     list(.) %>%
  #     append(read_rds("models/kfold/q2q3_m2_llk.rds"), .) %>%
  #     write_rds(., "models/kfold/q2q3_m2_llk.rds")
  #   }
  # }
  
  # # Model 3
  # for (k in 1:10) {
  #   print(paste("[", k, "/10] - Model 3", sep = ""))
  #   fit <- stan(
  #     "models/kfold/q2q3_m3_kfold.stan",
  #     data = data_fold(k), 
  #     seed = seeds[4],
  #     iter = 1000, 
  #     warmup = 750, 
  #     pars = c("log_lik_h")
  #   )
  #   if (k == 1) {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>% 
  #     list(.) %>%
  #     write_rds(., "models/kfold/q2q3_m3_llk.rds")
  #   } else {
  #     extract_log_lik(fit, parameter_name = "log_lik_h") %>%
  #     list(.) %>%
  #     append(read_rds("models/kfold/q2q3_m3_llk.rds"), .) %>%
  #     write_rds(., "models/kfold/q2q3_m3_llk.rds")
  #   }
  # }
  
  # Model 4
  for (k in 5:10) {
    print(paste("[", k, "/10] - Model 4", sep = ""))
    fit <- stan(
      "models/kfold/q2q3_m4_kfold.stan",
      data = data_fold(k), 
      seed = seeds[5],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h")
    )
    if (k == 1) {
      extract_log_lik(fit, parameter_name = "log_lik_h") %>% 
      list(.) %>%
      write_rds(., "models/kfold/q2q3_m4_llk.rds")
    } else {
      extract_log_lik(fit, parameter_name = "log_lik_h") %>%
      list(.) %>%
      append(read_rds("models/kfold/q2q3_m4_llk.rds"), .) %>%
      write_rds(., "models/kfold/q2q3_m4_llk.rds")
    }
  }
  
  # Model 5
  for (k in 1:10) {
    print(paste("[", k, "/10] - Model 5", sep = ""))
    fit <- stan(
      "models/kfold/q2q3_m5_kfold.stan",
      data = data_fold(k), 
      seed = seeds[6],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h")
    )
    if (k == 1) {
      extract_log_lik(fit, parameter_name = "log_lik_h") %>% 
      list(.) %>%
      write_rds(., "models/kfold/q2q3_m5_llk.rds")
    } else {
      extract_log_lik(fit, parameter_name = "log_lik_h") %>%
      list(.) %>%
      append(read_rds("models/kfold/q2q3_m5_llk.rds"), .) %>%
      write_rds(., "models/kfold/q2q3_m5_llk.rds")
    }
  }
  
  # Model 6
  for (k in 1:10) {
    print(paste("[", k, "/10] - Model 6", sep = ""))
    fit <- stan(
      "models/kfold/q2q3_m6_kfold.stan",
      data = data_fold(k), 
      seed = seeds[7],
      iter = 1000, 
      warmup = 750, 
      pars = c("log_lik_h")
    )
    if (k == 1) {
      extract_log_lik(fit, parameter_name = "log_lik_h") %>% 
      list(.) %>%
      write_rds(., "models/kfold/q2q3_m6_llk.rds")
    } else {
      extract_log_lik(fit, parameter_name = "log_lik_h") %>%
      list(.) %>%
      append(read_rds("models/kfold/q2q3_m6_llk.rds"), .) %>%
      write_rds(., "models/kfold/q2q3_m6_llk.rds")
    }
  }
  

# Estimate ----------------------------------------------------------------
  m0_llk <- read_rds("models/kfold/q2q3_m0_llk.rds")
  m1_llk <- read_rds("models/kfold/q2q3_m1_llk.rds")
  m2_llk <- read_rds("models/kfold/q2q3_m2_llk.rds")
  m3_llk <- read_rds("models/kfold/q2q3_m3_llk.rds")
  m4_llk <- read_rds("models/kfold/q2q3_m4_llk.rds")
  m5_llk <- read_rds("models/kfold/q2q3_m5_llk.rds")
  m6_llk <- read_rds("models/kfold/q2q3_m6_llk.rds")

  # Calculate expected log pointwise predictive density for each model
  llk    <- list(m0_llk, m1_llk, m2_llk, m3_llk, m4_llk, m5_llk, m6_llk)
  elpd_i <- map(llk, ~map(., ~log(colMeans(exp(.))))) %>%
            map(., unlist)
  
  # Compare models
  elpd_kfold <- tibble(
    model  = 0:(length(elpd_i)-1),
    elpd_i = elpd_i
  ) %>% 
    mutate(
      elpd    = map(elpd_i, sum) %>% unlist(),
      elpd_se = map(elpd_i, ~(sd(.) * sqrt(length(.)))) %>% unlist(),
      elpd_d  = map2(lag(elpd_i), elpd_i, compare_kfold)
    ) %>% unnest(elpd_d)
  
  # Export
  elpd_kfold %>% write_rds(., "models/kfold/q2q3_elpd_kfold.rds")
  
  # Calculate stacking weights
  matrix(unlist(elpd_i), ncol = length(elpd_i)) %>% stacking_weights()
  matrix(unlist(elpd_i), ncol = length(elpd_i)) %>% pseudobma_weights()
