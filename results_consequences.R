rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(tidybayes); library(loo)

  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seeds = c(1157744, 2691555, 4045395, 2910499, 5195618, 5975015, 2391421)
  

# Functions ---------------------------------------------------------------


# Import ------------------------------------------------------------------
  dl <- read_rds("data/d1.rds") %>%
        filter(!is.na(q1)) %>%
        mutate_at(vars(q2, q3), funs(mean, sd), na.rm = TRUE) %>%
        mutate(q2_z = (q2 - q2_mean) / q2_sd, q3_z = (q3 - q3_mean) / q3_sd) %>%
        arrange(participant, target)
  
  # Prepare for imputation
  dl_imp <- dl %>% 
            select(participant, target, q2_z, q3_z) %>%
            gather("item", "y", q2_z, q3_z) %>%
            mutate(ii = 1:n())

  # Compose data list
  data_list <- with(dl, list(
    # Numbers
    N_row = length(participant),
    N_col = 2,
    J = n_distinct(participant),
    K = n_distinct(category),
    N_obs = sum(!is.na(dl_imp$y)),
    N_mis = sum(is.na(dl_imp$y)),
    N_all_obs = dl %>% filter(!is.na(q2), !is.na(q3)) %>% nrow(),
    # Indices
    jj = participant,
    kk = category,
    ii_obs = dl_imp$ii[!is.na(dl_imp$y)],
    ii_mis = dl_imp$ii[is.na(!dl_imp$y)],
    nn_obs = dl %>% mutate(nn = 1:n()) %>% filter(!is.na(q2), !is.na(q3)) %>% pull(nn),
    # Vectors
    x_scst = ifelse(ig_caste == "scst", 1L, 0L),
    x_obc = ifelse(ig_caste == "obc", 1L, 0L),
    x_q1 = q1,
    y_obs = dl_imp$y[!is.na(dl_imp$y)]
  ))


# Models ------------------------------------------------------------------
  # Model: Varying intercept (Participant)
  m0_stan <- stan("models/q2q3_m0.stan", data = data_list, seed = seeds[1],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)

  # Model: M0 + Varying Intercept (Target: Category)
  m1_stan <- stan("models/q2q3_m1.stan", data = data_list, seed = seeds[2],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)

  # Model: M1 + Participant: SCST * Target: Category (4)
  m2_stan <- stan("models/q2q3_m2.stan", data = data_list, seed = seeds[3],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)

  # Model: M2 + Participant: OBC * Target: Category (4)
  m3_stan <- stan("models/q2q3_m3.stan", data = data_list, seed = seeds[4],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)
  
  # Model: M3 + Participant: Q1
  m4_stan <- stan("models/q2q3_m4.stan", data = data_list, seed = seeds[5],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)
  
  # Model: M3 + Participant: Q1  * Target: Category (6)
  m5_stan <- stan("models/q2q3_m5.stan", data = data_list, seed = seeds[6],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)
  
  # Model: M5 + Participant: Q1 * Participant: SCST/OBC * Target: Category (4)
  m6_stan <- stan("models/q2q3_m6.stan", data = data_list, seed = seeds[7],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)
  
  
# Compare -----------------------------------------------------------------
  # Compare models
  mc_fit <- list(m0_stan, m1_stan, m2_stan, m3_stan, m4_stan, m5_stan, m6_stan)
  mc_R2  <- map(mc_fit, ~gather_draws(., R2[item]) %>% median_hdi(.width = .97))

  # Export
  write_rds(mc_fit, "models/q2q3_mc_fit.rds")

  
# Predict -----------------------------------------------------------------
  # Posterior draws
  draws <- left_join(
      spread_draws(m5_stan, b_0[item], b_k[category,item], b_q1[category,item], 
                   sigma_k[item], sigma_k[item], sigma_scst_k[item], 
                   sigma_obc_k[item], sigma_q1_k[item]),
      spread_draws(m5_stan, b_scst_k[category,item], b_obc_k[category,item])
    ) %>% replace_na(
      replace = list(b_scst_k = 0, b_obc_k = 0)
    ) %>% select(
      -.chain, -.iteration
    )
  
  # Posterior predictions
  post <- crossing(
    category = 1:6,
    ig_caste = 1:3,
    item = 1:2,
    x_q1 = 0:1
  ) %>% mutate(
    x_obc  = ifelse(ig_caste == 2, 1, 0),
    x_scst = ifelse(ig_caste == 3, 1, 0)
  ) %>% full_join(
    draws, by = c("item", "category")
  ) %>% mutate(
    m_est_z = b_0 + b_k*sigma_k + b_q1*sigma_q1_k*x_q1 + b_scst_k*sigma_scst_k*
              x_scst + b_obc_k*sigma_obc_k*x_obc,
    item = recode(item, "1" = "q2", "2" = "q3")
  )
  
  # Unstandardize
  post <- dl %>% 
    distinct(q2_mean, q3_mean, q2_sd, q3_sd) %>% 
    gather("item", "response") %>% separate(item, c("item", "summary"), "_") %>% 
    spread(summary, response) %>%
    left_join(post, ., by = "item") %>%
    mutate(
      m_est = m_est_z * sd + mean
    ) %>%
    select(.draw, category:ig_caste, x_q1:x_scst, item, m_est, m_est_z)

  
  # Export
  write_rds(post, "models/q2q3_post.rds")
