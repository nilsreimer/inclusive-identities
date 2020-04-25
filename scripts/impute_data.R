rm(list = ls())

# Notes -------------------------------------------------------------------

  #########################################################################
  # This script imputes missing responses (n = 66) for the intergroup     #
  # contact variables, and returns a data frame with both observed        #
  # responses and draws from the posterior distribution for missing       #
  # responses.                                                            # 
  #########################################################################

# Library -----------------------------------------------------------------
  
  # Load packages
  library(tidyverse); library(rstan)

  # Stan options
  n_cores <- 8L
  options(mc.cores = n_cores)
  rstan_options(auto_write = TRUE)
  seed = 70475734


# Prepare -----------------------------------------------------------------

  # Import data
  d2 <- read_rds("data/d2.rds") %>% select(participant, ig_caste, q4:q23)
  
  # Prepare data (long format)
  dl <- d2 %>%
    pivot_longer(
      -participant:-ig_caste,
      names_to = "item",
      names_prefix = "q",
      names_ptypes = list(item = integer()),
      values_to = "x_obs"
    ) %>% 
    arrange(item, participant) %>%
    rowid_to_column("ii") %>% 
    select(ii, participant, ig_caste, item, x_obs)
  
  # Prepare data (wide format)
  dw <- d2 %>%
    mutate(
      kk = as.integer(factor(ig_caste))
    ) %>%
    select(participant, ig_caste, kk) %>%
    arrange(participant)
  
  # Standardize
  dl <- dl %>% 
    group_by(ig_caste, item) %>% 
    mutate(
      z_obs = (x_obs - mean(x_obs, na.rm = TRUE))/sd(x_obs, na.rm = TRUE)
    ) %>% 
    ungroup()
  
  # Compose data list
  data_list <- with(dl, list(
      N_row = n_distinct(participant),
      N_col = n_distinct(item),
      N_obs = sum(!is.na(x_obs)),
      N_mis = sum(is.na(x_obs)),
      ii_obs = ii[!is.na(x_obs)],
      ii_mis = ii[is.na(x_obs)],
      x_obs = z_obs[!is.na(x_obs)]
    )) %>% 
    append(
      with(dw, list(
        K = n_distinct(kk),
        kk = kk
      ))
    )
  

# Run ---------------------------------------------------------------------

  # Run model
  fit <- stan(
    "models/impute_data.stan",
    data = data_list,
    iter = 1000 + 4000/n_cores, 
    warmup = 1000,
    chains = n_cores,
    init = function(chain_id) list(
      x_mu = array(0, dim = c(3, 20)),
      x_sigma = array(1, dim = c(3, 20))
    ),
    seed = seed
  )
  
  # Inspect
  check_hmc_diagnostics(fit)
  stan_rhat(fit, pars = c("x_mu", "x_sigma", "L_Omega", "x_mis"))
  print(fit, pars = "x_mu")
  print(fit, pars = "x_sigma")
  print(fit, pars = "Rho")
  print(fit, pars = "x_mis")


# Extract -----------------------------------------------------------------

  # Extract posterior draws
  dl_imp <- as.data.frame(fit, pars = "x_imp") %>% 
    rowid_to_column(".draw") %>% 
    pivot_longer(
      -.draw,
      names_to = "ii",
      names_pattern = "x_imp\\[([0-9]*)\\]",
      names_ptypes = list(ii = integer()),
      values_to = "z_imp"
    )
  
  # Merge with observed responses
  dl_imp <- dl_imp %>% 
    left_join(dl, by = "ii") %>% 
    group_by(.draw, ig_caste, item) %>% 
    mutate(
      x_imp = z_imp * sd(x_obs, na.rm = TRUE) + mean(x_obs, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    mutate(
      category = case_when(
        item %in% 4:8 ~ 3L,
        item %in% 9:13 ~ 2L,
        item %in% 14:18 ~ 1L,
        item %in% 19:23 ~ 4L
      ),
      item = case_when(
        item %in% c(4L,  9L, 14L, 19L) ~ "cq",
        item %in% c(5L, 10L, 15L, 20L) ~ "pc",
        item %in% c(6L, 11L, 16L, 21L) ~ "nc",
        item %in% c(7L, 12L, 17L, 22L) ~ "of",
        item %in% c(8L, 13L, 18L, 23L) ~ "of",
      ),
      item = factor(item, levels = c("cq", "pc", "nc", "of"))
    ) %>% 
    select(participant, ig_caste, category, item, x_obs, z_obs, .draw, x_imp, z_imp) %>% 
    arrange(.draw, participant, category, item)
  
  # Simplify
  dl_imp <- dl_imp %>% 
    mutate(.draw = if_else(is.na(x_obs), .draw, NA_integer_)) %>% 
    distinct()

  # Export
  write_rds(dl_imp, "results/ic_impute.rds")
