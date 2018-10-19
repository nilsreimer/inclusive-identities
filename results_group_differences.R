rm(list = ls())

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes); library(rstan); library(loo)

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


# Models ------------------------------------------------------------------
  # Model: Varying Intercept (Participant)
  m0_stan <- stan("models/q1_m0.stan", data = data_list, seed = seeds[1],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Model: + Varying Intercept (Target: Category)
  m1_stan <- stan("models/q1_m1.stan", data = data_list, seed = seeds[2],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Model: + Participant: SCST * Target: Category (4)
  m2_stan <- stan("models/q1_m2.stan", data = data_list, seed = seeds[3],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Model: + Participant: OBC * Target: Category (4)
  m3_stan <- stan("models/q1_m3.stan", data = data_list, seed = seeds[4],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)


# Compare -----------------------------------------------------------------
  # Compare models
  mc_fit <- list(m0_stan, m1_stan, m2_stan, m3_stan)
  
  # Export
  write_rds(mc_fit, "models/q1_mc_fit.rds")


# Predict -----------------------------------------------------------------
  # Extract posterior samples
  post <- tibble(
    model = 0:3,
    fit = mc_fit
  ) %>% mutate(
    p_est = map(fit, ~spread_draws(., p_k[category, ig_caste]))
  ) %>% select(model, p_est)
  
  # Export
  write_rds(post, "models/q1_post.rds")
