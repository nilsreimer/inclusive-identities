rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------

  # Load packages
  library(tidyverse); library(rstan); library(loo); library(tidybayes)
  
  # Stan options
  n_cores <- 8L
  options(mc.cores = n_cores)
  rstan_options(auto_write = TRUE)
  seeds = c(9658217, 4125367, 9337593, 1407167)

  
# Prepare -----------------------------------------------------------------

  # Import data
  dw <- read_rds("data/d2.rds") %>%
    transmute(
      participant,
      ig_caste,
      ps_scst     = (q25 + q26)/2,
      ps_obc      = (q27 + q28)/2,
      ps_muslim   = q29,
      x_scst      = if_else(ig_caste == "scst", 1L, 0L),
      x_obc       = if_else(ig_caste == "obc", 1L, 0L),
      x_gm        = if_else(ig_caste == "gm", 1L, 0L),
      x_nc_gm     = q16,
      x_of_gm     = (q17 + q18)/2,
      x_nc_obc    = q11,
      x_of_obc    = (q12 + q13)/2,
      x_nc_scst   = q6,
      x_of_scst   = (q7 + q8)/2,
      x_nc_muslim = q21,
      x_of_muslim = (q22 + q23)/2,
    ) %>%
    left_join(
      read_rds("data/d1.rds") %>% 
        select(participant, category, q1) %>% 
        filter(category <= 4L) %>% 
        mutate(category = recode(
          category, "1" = "gm", "2" = "obc", "3" = "scst", "4" = "muslim"
        )) %>% 
        group_by(participant, category) %>% 
        summarise(q1 = mean(q1)) %>% 
        ungroup() %>% 
        pivot_wider(
          names_from = category,
          names_prefix = "x_q1_",
          values_from = q1
        ),
      by = "participant"
    ) %>% 
    mutate_at(vars(starts_with("x_")), as.double) %>% 
    mutate_at(
      vars(x_nc_gm, x_of_gm, x_q1_gm), 
      ~if_else(ig_caste == "gm", NA_real_, .)
    ) %>% 
    mutate_at(
      vars(x_nc_obc, x_of_obc, x_q1_obc), 
      ~if_else(ig_caste == "obc", NA_real_, .)
    ) %>% 
    mutate_at(
      vars(x_nc_scst, x_of_scst, x_q1_scst),
      ~if_else(ig_caste == "scst", NA_real_, .)
    ) %>% 
    arrange(participant)
  
  # Standardize
  dz <- dw %>% 
    mutate_at(
      vars(starts_with("ps_"), starts_with("x_nc_"), starts_with("x_of_")), 
      ~(. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)
    ) %>% 
    mutate_at(
      vars(starts_with("x_q1_")),
      ~(. - mean(., na.rm = TRUE))
    )
  
  # Exclude missing observations
  dz <- dz %>% 
    mutate_at(
      vars(x_nc_gm, x_of_gm, x_q1_gm),
      ~if_else(ig_caste == "gm", -99, .)
    ) %>% 
    mutate_at(
      vars(x_nc_obc, x_of_obc, x_q1_obc),
      ~if_else(ig_caste == "obc", -99, .)
    ) %>% 
    mutate_at(
      vars(x_nc_scst, x_of_scst, x_q1_scst),
      ~if_else(ig_caste == "scst", -99, .)
    ) %>% 
    drop_na()

  # Compose data list
  data_list <- with(dz, list(
    N = length(participant),
    K = 3L,
    y = cbind(ps_scst, ps_obc, ps_muslim)
  ))
  
  
# Run ---------------------------------------------------------------------
  
  # Model: Group differences (categories)
  m0_fit <- stan(
      "models/ps_m0.stan", 
      data = data_list, 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      seed = seeds[1],
      pars = c("mu"), include = FALSE
    )
  
  # Model: + Fixed effects (SC/ST + OBC)
  m1_fit <- data_list %>% 
    append(with(dz, list(
      x_obc  = x_obc,
      x_scst = x_scst
    ))) %>% 
    stan(
      "models/ps_m1.stan", 
      data = ., 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      seed = seeds[2],
      pars = c("mu"), include = FALSE
    )
  
  # Model: + Fixed effects (NC + OF)
  m2_fit <- data_list %>% 
    append(with(dz, list(
      x_gm = x_gm, 
      x_obc = x_obc,
      x_scst = x_scst,
      x_nc_gm = x_nc_gm,
      x_of_gm = x_of_gm,
      x_nc_obc = x_nc_obc,
      x_of_obc = x_of_obc,
      x_nc_scst = x_nc_scst,
      x_of_scst = x_of_scst,
      x_nc_muslim = x_nc_muslim,
      x_of_muslim = x_of_muslim
    ))) %>% 
    stan(
      "models/ps_m2.stan", 
      data = ., 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      seed = seeds[3],
      pars = c("mu"), include = FALSE
    )
  
  # Model: + Fixed effects (categorization)
  m3_fit <- data_list %>% 
    append(with(dz, list(
      x_gm = x_gm,
      x_obc = x_obc,
      x_scst = x_scst,
      x_q1_gm = x_q1_gm,
      x_q1_obc = x_q1_obc,
      x_q1_scst = x_q1_scst,
      x_q1_muslim = x_q1_muslim
    ))) %>% 
    stan(
      "models/ps_m3.stan", 
      data = ., 
      iter = 1000 + 8000/n_cores, warmup = 1000, 
      chains = n_cores,
      seed = seeds[4],
      pars = c("mu"), include = FALSE
    )
  

# Compare -----------------------------------------------------------------

  # Compile model results
  ps_fit <- tibble(
      model = 0:3,
      fit = list(m0_fit, m1_fit, m2_fit, m3_fit)
    ) %>% 
    mutate(
      R2 = map_dbl(fit, ~as.data.frame(., pars = "R2") %>% rowMeans() %>% median()),
      log_lik = map(fit, ~extract_log_lik(., merge_chains = FALSE)),
      elpd = map(log_lik, ~loo(., r_eff = relative_eff(exp(.)))),
      w = loo_model_weights(elpd, method = "pseudobma")
    )


# Extract -----------------------------------------------------------------

  # Extract coefficients
  post <- m1_fit %>% 
    spread_draws(b_0[kk], b_obc[kk], b_scst[kk]) %>% 
    ungroup() %>% 
    crossing(ig_caste = c("gm", "obc", "scst")) %>% 
    mutate(
      m_est = case_when(
        ig_caste == "gm" ~ b_0,
        ig_caste == "obc" ~ b_0 + b_obc,
        ig_caste == "scst" ~ b_0 + b_scst
      ),
      category = recode(
        kk,
        "1" = "scst",
        "2" = "obc",
        "3" = "muslim"
      )
    ) %>% 
    left_join(
      dw %>% 
        select(starts_with("ps_")) %>% 
        pivot_longer(
          starts_with("ps_"),
          names_to = "category",
          names_prefix = "ps_",
          values_to = "value"
        ) %>% 
        group_by(category) %>% 
        summarize_at(vars(value), list("mean" = mean, "sd" = sd), na.rm = TRUE),
      by = "category"
    ) %>% 
    mutate(m_est = m_est * sd + mean) %>% 
    select(.chain:.draw, kk, ig_caste, category, m_est) %>% 
    arrange(.draw, kk, ig_caste)
  
  # Export
  write_rds(post, "results/ps_post.rds")
  