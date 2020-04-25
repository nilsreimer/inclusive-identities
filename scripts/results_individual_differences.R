rm(list = ls())

# Notes -------------------------------------------------------------------

  #########################################################################
  # Adjust n_cores to the number of CPU cores you want to use.            #
  #########################################################################

# Library -----------------------------------------------------------------

  # Load packages
  library(tidyverse); library(tidybayes); library(rstan)
  
  # Stan options
  n_cores <- 8L
  options(mc.cores = n_cores)
  rstan_options(auto_write = TRUE)
  seeds = c(6678002, 8970252,	9556204, 8660009)
  

# Functions ---------------------------------------------------------------
  inv_logit <- function(x) exp(x) / ( 1 + exp(x) )
  
  
# Prepare -----------------------------------------------------------------

  # Import data
  dl <- read_rds("data/d1.rds") %>%
    filter(!is.na(q1)) %>%
    mutate(og = case_when(
      category == 1 & ig_caste == "gm" ~ 0L,
      category == 2 & ig_caste == "obc" ~ 0L,
      category == 3 & ig_caste == "scst" ~ 0L,
      TRUE ~ 1L
    )) %>%
    arrange(participant, target)
  
  # Import imputed data
  dl_imp <- read_rds("results/ic_impute.rds") %>% 
    mutate(og = case_when(
      category == 1 & ig_caste == "gm" ~ 0L,
      category == 2 & ig_caste == "obc" ~ 0L,
      category == 3 & ig_caste == "scst" ~ 0L,
      TRUE ~ 1L
    )) %>% 
    select(-z_obs, -z_imp)
  
  # Summarize/Standardize
  dl_imp <- dl_imp %>% 
    filter(!is.na(x_obs)) %>% 
    select(-.draw) %>% 
    crossing(.draw = 1:max(dl_imp$.draw, na.rm = TRUE)) %>% 
    bind_rows(dl_imp %>% filter(is.na(x_obs))) %>% 
    group_by_at(vars(-x_obs, -x_imp)) %>% 
    summarise_at(vars(x_obs, x_imp), mean) %>% 
    group_by(.draw, item) %>% 
    mutate(
      x_obs = if_else(og == 0L, NA_real_, x_obs),
      z_obs = (x_obs - mean(x_obs, na.rm = TRUE))/sd(x_obs, na.rm = TRUE),
      z_imp = (x_imp - mean(x_obs, na.rm = TRUE))/sd(x_obs, na.rm = TRUE)
    ) %>% 
    ungroup()
  
  # Summarize posterior distributions for model
  dl_imp <- dl_imp %>% 
    group_by(participant, category, item) %>% 
    summarize(
      x_obs = unique(x_obs),
      z_obs = unique(z_obs),
      z_mu = mean(z_imp),
      z_sigma = sd(z_imp)
    ) %>% 
    ungroup()
  
  # Merge missing and observed data
  dl_imp <- dl_imp %>% 
    complete(participant, item, category = 1:6) %>% 
    left_join(
      dl %>% select(participant, target, category, og), 
      .,
      by = c("participant", "category")
    ) %>% 
    mutate(
      z_obs = if_else(og == 0L | category > 4L, -99, z_obs)
    )
  
  # Arrange for model
  dl_imp <- dl_imp %>% 
    arrange(item, participant, target) %>% 
    rowid_to_column("ii")

  # Add factor scores for SDO-D and SDO-E
  source("scripts/prepare_data_sdo.R")

  # Compose data list
  data_list <- with(dl, list(
      N = length(participant),
      J = n_distinct(participant),
      K = n_distinct(category),
      jj = participant,
      kk = category,
      x_scst = ifelse(ig_caste == "scst", 1L, 0L),
      x_obc = ifelse(ig_caste == "obc", 1L, 0L),
      x_og = og,
      x_gs = gs,
      x_sdo_d = sdo_d,
      x_sdo_e = sdo_e,
      y = q1
    )) %>% 
    append(
      with(dl_imp, list(
        N_obs = sum(!is.na(z_obs)),
        N_mis = sum(is.na(z_obs)),
        ii_obs = ii[!is.na(z_obs)],
        ii_mis = ii[is.na(z_obs)],
        x_obs = z_obs[!is.na(z_obs)],
        x_mu = z_mu[is.na(z_obs)],
        x_sigma = z_sigma[is.na(z_obs)]
      ))
    )


# Run ---------------------------------------------------------------------

  # Model: M2 + Participant: Intergroup contact (4)
  m4_fit <- stan(
    "models/q1_m4.stan", 
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[1],
    pars = c("alpha", "x_imp", "x_cq", "x_pc", "x_nc", "x_of"), 
    include = FALSE
  )
  
  # Inspect
  m4_fit %>% 
    gather_draws(b_cq, b_pc, b_nc, b_of) %>% 
    mutate(.value = exp(.value)) %>% 
    median_qi(.value, .width = 0.97) %>% 
    mutate_if(is.double, round, digits = 2)
  
  # Model: M2 + Participant: Intergroup contact (2)
  m5_fit <- stan(
    "models/q1_m5.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    seed = seeds[2],
    pars = c("alpha", "x_imp", "x_nc", "x_of"),
    include = FALSE
  )

  # Inspect
  m5_fit %>%
    gather_draws(b_nc, b_of) %>%
    mutate(.value = exp(.value)) %>% 
    median_qi(.value, .width = 0.97) %>% 
    mutate_if(is.double, round, digits = 2)
  
  # Model: M5 + Participant: Intergroup contact (2) * Target: Category (4)
  m6_fit <- stan(
    "models/q1_m6.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[3],
    pars = c("alpha", "x_imp", "x_nc", "x_of"), 
    include = FALSE
  )
  
  # Inspect
  m6_fit %>%
    spread_draws(b_nc, b_nc_k[kk], b_of, b_of_k[kk]) %>%
    mutate(b_nc = b_nc + b_nc_k, b_of = b_of + b_of_k) %>%
    mutate_if(is.double, exp) %>%
    median_qi(b_nc, b_pc)

  # Model: M2 + Participant: SDO (2)
  m7_fit <- stan(
    "models/q1_m7.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[4],
    pars = c("alpha", "x_imp", "x_nc", "x_of"),
    include = FALSE
  )
  
  # Inspect
  m7_fit %>% 
    gather_draws(b_sdo_d, b_sdo_e) %>% 
    mutate(.value = exp(.value)) %>% 
    median_qi(.value, .width = 0.97) %>% 
    mutate_if(is.double, round, digits = 2)


# Extract -----------------------------------------------------------------

  # Compile model results
  q1_fit <- list(m4_fit, m5_fit, m_6_fit, m7_fit)

  # Merge with m0_stan to m3_stan
  q1_fit <- append(read_rds("results/q1_fit.rds"), q1_fit)

  # Export
  write_rds(q1_fit, "results/q1_fit.rds")
  
  # Extract posterior samples
  post <- dl_imp %>% 
    filter(og == 1L, category <= 4, item %in% c("nc", "of")) %>% 
    distinct(participant, category, item, x_obs) %>% 
    group_by(item) %>% 
    summarise(
      mean = mean(x_obs, na.rm = TRUE),
      sd = sd(x_obs, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    pivot_wider(
      names_from = item,
      values_from = c(mean, sd)
    ) %>% 
    crossing(
      category = 1:4,
      ig_caste = 1:3,
      nesting(
        item = c(rep("nc", 5), rep("of", 9)),
        response = c(seq(1, 5, 1), seq(1, 5, 0.5))
      )
    ) %>% 
    mutate(
      x_nc = if_else(
        item == "nc",
        (response - mean_nc)/sd_nc,
        0
      ),
      x_of = if_else(
        item == "of",
        (response - mean_of)/sd_of,
        0
      )
    ) %>% 
    select(-mean_nc:-sd_of) %>% 
    filter(category != ig_caste) %>% 
    full_join(
      m5_fit %>%
        spread_draws(l_k[category, ig_caste], b_nc, b_of) %>%
        filter(category <= 4, category != ig_caste),
      by = c("category", "ig_caste")
    ) %>%
    mutate(
      p_est = inv_logit(l_k + b_nc * x_nc + b_of * x_of)
    ) %>%
    select(.draw, category, ig_caste, item, response, x_nc, x_of, p_est) %>%
    ungroup()
  
  # Export
  write_rds(post, "results/q1_ic_post.rds")
