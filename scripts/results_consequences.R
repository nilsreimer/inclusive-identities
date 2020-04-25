rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------
  
  # Load packages
  library(tidyverse); library(rstan); library(tidybayes)

  # Stan options
  n_cores <- 8L
  options(mc.cores = n_cores)
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
  
  # Prepare for imputation
  dl_imp <- dl %>% 
    select(participant, target, q2_z, q3_z) %>%
    gather("item", "y", q2_z, q3_z) %>%
    arrange(item, participant, target) %>% 
    rowid_to_column("ii")
    
  # Compose data list
  data_list <- with(dl, list(
      N_row = length(participant),
      N_col = 2,
      J = n_distinct(participant),
      K = n_distinct(category),
      jj = participant,
      kk = category,
      x_scst = if_else(ig_caste == "scst", 1L, 0L),
      x_obc = if_else(ig_caste == "obc", 1L, 0L),
      x_q1 = q1
    )) %>% 
    append(
      with(dl_imp, list(
        N_obs = sum(!is.na(y)),
        N_mis = sum(is.na(y)),
        ii_obs = ii[!is.na(y)],
        ii_mis = ii[is.na(y)],
        y_obs = y[!is.na(y)]
      ))
    )


# Models ------------------------------------------------------------------
  
  # Model: Varying intercept (Participant)
  m0_fit <- stan(
    "models/q2q3_m0.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    seed = seeds[1],
    pars = c("Mu", "Y", "y", "b_j_z"),
    include = FALSE
  )

  # Inspect
  print(m0_fit, pars = c("b_0", "sigma", "Rho", "sigma_j", "Rho_j", "R2"))

  # Model: M0 + Varying Intercept (Target: Category)
  m1_fit <- stan(
    "models/q2q3_m1.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[2],
    pars = c("Mu", "Y", "y", "b_j_z"),
    include = FALSE
  )

  # Inspect
  print(m1_fit, pars = c(
    "b_0", "sigma", "Rho",
    "sigma_j", "Rho_j", "sigma_k",
    "R2"
  ))

  # Model: M1 + Participant: SCST * Target: Category (4)
  m2_fit <- stan(
    "models/q2q3_m2.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[3],
    pars = c("Mu", "Y", "y", "b_j_z"),
    include = FALSE
  )

  # Inspect
  print(m2_fit, pars = c(
    "b_0", "b_scst", "sigma", "Rho", "sigma_j", "Rho_j",
    "sigma_k", "sigma_scst_k", "R2", "b_scst_k_z"
  ))

  # Model: M2 + Participant: OBC * Target: Category (4)
  m3_fit <- stan(
    "models/q2q3_m3.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[4],
    pars = c("Mu", "Y", "y", "b_j_z"),
    include = FALSE
  )

  # Inspect
  print(m3_fit, pars = c(
    "b_0", "b_scst", "b_obc", "sigma", "Rho", "sigma_j", "Rho_j",
    "sigma_k", "sigma_scst_k", "sigma_obc_k", "R2", "b_obc_k_z"
  ))

  # Model: M3 + Fixed effect (Categorization)
  m4_fit <- stan(
    "models/q2q3_m4.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[5],
    pars = c("Mu", "Y", "y", "b_j_z"),
    include = FALSE
  )

  # Inspect
  print(m4_fit, pars = c(
    "b_0", "b_scst", "b_obc", "b_q1", "sigma", "Rho", "sigma_j", "Rho_j",
    "sigma_k", "sigma_scst_k", "sigma_obc_k", "R2", "b_obc_k_z"
  ))

  # Model: M4 + Varying effect (Categories)
  m5_fit <- stan(
    "models/q2q3_m5.stan", 
    data = data_list, 
    iter = 1000 + 8000/n_cores, warmup = 1000, 
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[6],
    pars = c("Mu", "Y", "y", "b_j_z"), 
    include = FALSE
  )
  
  # Inspect
  print(m5_fit, pars = c(
    "b_0", "b_scst", "b_obc", "sigma", "Rho", "sigma_j", "Rho_j", "R2",
    "sigma_k", "sigma_scst_k", "sigma_obc_k", "b_obc_k_z"
  ))
  
  # Model: M5 + Interaction (SC/ST, OBC)
  m6_fit <- stan(
    "models/q2q3_m6.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[7],
    pars = c("Mu", "Y", "y", "b_j_z"),
    include = FALSE
  )

  # Inspect
  print(m6_fit, pars = c(
    "b_0", "b_scst", "b_obc", "b_q1", "b_q1_k_z", "sigma_q1_k",
    "b_q1_scst", "b_q1_scst_k_z", "sigma_q1_scst_k",
    "b_q1_obc", "b_q1_obc_k_z", "sigma_q1_obc_k",
    "sigma", "Rho", "sigma_j", "Rho_j", "R2",
    "sigma_k", "sigma_scst_k", "sigma_obc_k", "b_obc_k_z"
  ))


# Extract -----------------------------------------------------------------
  
  # Compile model results
  q2q3_fit <- list(m0_fit, m1_fit, m2_fit, m3_fit, m4_fit, m5_fit, m6_fit)

  # Export
  write_rds(q2q3_fit, "results/q2q3_fit.rds")
  
  # Extract posterior samples
  post <- m5_fit %>% 
    spread_draws(
      b_0[m], b_k_z[kk, m], sigma_k[m], 
      b_q1[m], b_q1_k_z[kk, m], sigma_q1_k[m]
    ) %>% 
    left_join(
      spread_draws(
        model = m5_fit,
        b_scst[m], b_scst_k_z[kk, m], sigma_scst_k[m],
        b_obc[m], b_obc_k_z[kk, m], sigma_obc_k[m]
      ),
      by = c(".chain", ".iteration", ".draw", "m", "kk")
    ) %>% 
    ungroup() %>% 
    crossing(
      x_q1 = 0:1,
      nesting(
        x_obc  = c(0L, 1L, 0L),
        x_scst = c(0L, 0L, 1L)
      )
    ) %>% 
    mutate(
      z_est = if_else(
        kk <= 4L,
        b_0 + b_k_z*sigma_k + x_q1*(b_q1 + b_q1_k_z*sigma_q1_k) + x_scst*(b_scst + b_scst_k_z*sigma_scst_k) + x_obc*(b_obc + b_obc_k_z*sigma_obc_k),
        b_0 + b_k_z*sigma_k + x_q1*(b_q1 + b_q1_k_z*sigma_q1_k)
      )
    )
  
  # Recode
  post <- post %>% 
    transmute(
      .draw,
      outcome = if_else(m == 1, "q2", "q3"),
      category = kk,
      ig_caste = case_when(
        x_obc == 1L ~ "OBC",
        x_scst == 1L ~ "SC/ST",
        TRUE ~ "GM"
      ),
      q1 = x_q1,
      z_est
    ) %>% 
    arrange(.draw, outcome, category, ig_caste, q1)

  # Unstandardize
  post <- dl %>% 
    select(ends_with("mean"), ends_with("sd")) %>% 
    distinct() %>% 
    pivot_longer(
      everything(),
      names_to = c("outcome", "summary"),
      names_sep = "_",
      values_to = "value"
    ) %>% 
    pivot_wider(names_from = summary, values_from = value) %>% 
    left_join(post, ., by = "outcome") %>% 
    mutate(
      m_est = z_est * sd + mean
    ) %>% 
    select(-mean, -sd)
  
  # Export
  write_rds(post, "results/q2q3_q1_post.rds")
