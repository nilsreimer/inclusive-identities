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
  seeds <- c(5477649, 3984665, 8989566, 9363790)


# Prepare -----------------------------------------------------------------

  # Import data
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
    x_scst = if_else(ig_caste == "scst", 1L, 0L),
    x_obc = if_else(ig_caste == "obc", 1L, 0L),
    y = q1
  ))


# Run ---------------------------------------------------------------------

  # Model: Varying Intercept (Participant)
  m0_fit <- stan(
    "models/q1_m0.stan", 
    data = data_list, 
    iter = 1000 + 8000/n_cores, warmup = 1000, 
    chains = n_cores,
    seed = seeds[1],
    pars = c("alpha"), include = FALSE
  )
  
  # Model: + Varying Intercept (Target: Category)
  m1_fit <- stan(
    "models/q1_m1.stan", 
    data = data_list, 
    iter = 1000 + 8000/n_cores, warmup = 1000, 
    chains = n_cores,
    control = list(adapt_delta = 0.90),
    seed = seeds[2],
    pars = c("alpha"), include = FALSE
  )
  
  # Model: + Participant: SCST * Target: Category (4)
  m2_fit <- stan(
    "models/q1_m2.stan",
    data = data_list,
    iter = 1000 + 8000/n_cores, warmup = 1000,
    chains = n_cores,
    seed = seeds[3],
    pars = c("alpha"), include = FALSE
  )
  
  # Model: + Participant: OBC * Target: Category (4)
  m3_fit <- stan(
    "models/q1_m3.stan", 
    data = data_list, 
    iter = 1000 + 8000/n_cores, warmup = 1000, 
    chains = n_cores,
    control = list(adapt_delta = 0.95),
    seed = seeds[4],
    pars = c("alpha"), include = FALSE
  )


# Extract -----------------------------------------------------------------

  # Compile model results
  q1_fit <- list(m0_fit, m1_fit, m2_fit, m3_fit)
  
  # Export
  write_rds(q1_fit, "results/q1_fit.rds")

  # Extract posterior samples
  post <- tibble(
      model = 0:3,
      fit = q1_fit
    ) %>% mutate(
      p_est = map(fit, ~spread_draws(., p_k[category, ig_caste]))
    ) %>% 
    select(model, p_est)

  # Export
  write_rds(post, "results/q1_post.rds")
