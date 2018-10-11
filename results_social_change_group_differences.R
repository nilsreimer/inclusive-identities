rm(list = ls())

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes); library(rstan); library(loo)

  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seeds = c(6096624, 7239687, 8749870, 6124014, 6071342,
            6890159, 9814684, 1014316, 5004873, 6578813)


# Import ------------------------------------------------------------------
  dw <- read_rds("data/d2.rds") %>%
    transmute(
      participant,
      ig_caste,
      ld_ig     = q24_1,
      ld_scst   = (q24_2 + q24_3)/2,
      ld_obc    = q24_4,
      ld_gm     = q24_5,
      ld_hindu  = q24_6,
      ld_muslim = q24_7,
      ps_scst   = (q25 + q26)/2,
      ps_obc    = (q27 + q28)/2,
      ps_muslim = q29
    ) %>%
    arrange(participant)
  
  # Prepare for imputation
  dl_imp <- dw %>% 
    gather("item", "y", -participant:-ig_caste) %>%
    mutate(ii = 1:n()) %>%
    group_by(item) %>%
    mutate_at(vars(y), funs(mean, sd), na.rm = TRUE) %>%
    mutate(y_z = (y - mean) / sd) %>%
    ungroup() %>%
    group_by(participant) %>% 
    mutate(missing = 9 - sum(!is.na(y))) %>%
    ungroup()
  
  # Compose data list
  data_list <- list(
    # Numbers
    N_row = length(dw$participant),
    N_col = 9,
    N_obs = sum(!is.na(dl_imp$y)),
    N_mis = sum(is.na(dl_imp$y)),
    N_all_obs = dl_imp %>% filter(missing == 0) %>% distinct(participant) %>% nrow(),
    # Indices
    ii_obs = dl_imp$ii[!is.na(dl_imp$y)],
    ii_mis = dl_imp$ii[is.na(!dl_imp$y)],
    nn_obs = dl_imp %>% filter(missing == 0) %>% distinct(participant) %>% pull(participant),
    # Vectors
    x_scst = ifelse(dw$ig_caste == "scst", 1L, 0L),
    x_obc = ifelse(dw$ig_caste == "obc", 1L, 0L),
    y_obs = dl_imp$y_z[!is.na(dl_imp$y)]
  )


# Models ------------------------------------------------------------------
  # Model: Intercept
  m0_stan <- stan("models/ldps_m0.stan", data = data_list, seed = seeds[1],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)
  
  # Model: Group differences
  m1_stan <- stan("models/ldps_m1.stan", data = data_list, seed = seeds[2],
                  iter = 2000, warmup = 1000,
                  pars = c("Mu", "Y", "y", "Lcorr"), include = FALSE)
  
# Compare -----------------------------------------------------------------
  # Compare models
  mc_fit <- list(m0_stan, m1_stan)
  mc_llk <- map(mc_fit, extract_log_lik, merge_chains = FALSE)
  mc_loo <- map(mc_llk, ~loo(., r_eff = relative_eff(exp(.)), cores = 4))
  
  # Compare
  compare(x = mc_loo)


# Predict -----------------------------------------------------------------
  # Extract posterior draws
  post <- m1_stan %>%
    spread_draws(
      b_0[item], b_scst[item], b_obc[item]
    ) %>% 
    ungroup() %>%
    transmute(
      .chain,
      .iteration,
      .draw,
      item,
      gm   = b_0,
      obc  = b_0 + b_obc,
      scst = b_0 + b_scst
    ) %>%
    gather(
      "ig_caste", "m_est", gm:scst
    ) 
  
  # Prepare for figures
  post <- post %>% 
    mutate(
      item = recode(
        item,
        "1" = "ld_ig",
        "2" = "ld_scst",
        "3" = "ld_obc",
        "4" = "ld_gm",
        "5" = "ld_hindu",
        "6" = "ld_muslim",
        "7" = "ps_scst",
        "8" = "ps_obc",
        "9" = "ps_muslim"
      )
    ) %>%
    left_join(
      dl_imp %>% distinct(item, mean, sd),
      by = "item"
    ) %>% 
    mutate(
      m_est = m_est * sd + mean
    )
  
  # Export
  write_rds(post, "models/ldps_post.rds")
    