rm(list = ls())

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes); library(rstan); library(loo)

  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seed = 0892561437


# Import ------------------------------------------------------------------
  dl <- read_rds("data/d1.rds") %>%
        filter(!is.na(q1))
  
  # Compose data list
  data_list <- with(dl, list(
    N = length(participant),
    J = n_distinct(participant),
    K = n_distinct(category),
    jj = participant,
    kk = category,
    scst = ifelse(ig_caste == "scst", 1L, 0L),
    obc = ifelse(ig_caste == "obc", 1L, 0L),
    y = q1
  ))


# Models ------------------------------------------------------------------
  # Model: Varying Intercept (Participant)
  m0_stan <- stan("models/q1_m0.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000, 
                  pars = c("alpha"), include = FALSE)
  
  # Model: + Varying Intercept (Target: Category)
  m1_stan <- stan("models/q1_m1.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000, 
                  pars = c("alpha"), include = FALSE)
  
  # Model: + Participant: SCST * Target: Category (4)
  m2_stan <- stan("models/q1_m2.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000, 
                  pars = c("alpha"), include = FALSE)
  
  # Model: + Participant: OBC * Target: Category (4)
  m3_stan <- stan("models/q1_m3.stan", data = data_list, seed = seed,
                  iter = 2000, warmup = 1000, 
                  pars = c("alpha"), include = FALSE)

  # Export
  mc_fit <- list(m0_stan, m1_stan, m2_stan, m3_stan)
  write_rds(mc_fit, "models/mc_fit_group_differences.rds")

# Compare -----------------------------------------------------------------
  mc_fit <- read_rds("models/mc_fit_group_differences.rds")
  mc_llk <- map(mc_fit, extract_log_lik, merge_chains = FALSE)
  mc_loo <- map(mc_llk, ~loo(., r_eff = relative_eff(exp(.)), cores = 4))

  compare(x = mc_loo)
  

# Predict -----------------------------------------------------------------
  post <- mc_fit[[3]] %>%
    spread_draws(p_k[category, ig_caste]) %>%
    median_hdi(.width = c(.67, .89, .97)) %>%
    mutate(ig_caste = recode(ig_caste, "1" = "gm", "2" = "obc", "3" = "scst"))
    
  pred <- mc_fit[[4]] %>% 
    spread_draws(p_k[category, ig_caste]) %>%
    median_hdi(.width = .97)

  pred %>%
    filter(.width == .97) %>% 
  ggplot(., aes(x = category, y = p_k, colour = ordered(ig_caste))) +
    geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = ordered(ig_caste)), colour = NA, alpha = 0.2) +
    scale_x_reverse(breaks = 1:6, minor_breaks = NULL) +
    scale_colour_viridis_d(begin = 0, end = 0.8) +
    scale_fill_viridis_d(begin = 0, end = 0.8) +
    coord_flip(ylim = c(0, 1)) +
    # facet_grid(. ~ ig_caste) +
    NULL
  