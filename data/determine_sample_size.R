rm(list = ls())

# Notes -------------------------------------------------------------------
  # Reanalysing data from van Dommelen et al.'s (2015) studies, we found
  # that 100 responses allowed reasonably precise estimates of the proportion
  # of "us" categorizations in each of eight target categories. Posterior
  # probabilities for these estimates had a precision of SD < 0.50 on the 
  # log odds scale.

  # This R script is included for review purposes only as we are not able to
  # make the underlying data available.


# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(tidybayes)
  
  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seeds = c(35702546)
  

# Import ------------------------------------------------------------------
  # Import wide format
  dw <- haven::read_spss("data/old/d_old.sav") %>%
        haven::zap_formats() %>%
        select(participant, card1:card24)
  
  # Spread to long format
  dl <- dw %>% 
        gather("target", "q1", card1:card24) %>%
        mutate(target = str_extract(target, "[0-9]+")) %>%
        mutate_at(vars(target, q1), as.integer) %>%
        arrange(participant, target)
  
  # Exclude completely missing cases
  dl <- dl %>%
        group_by(participant) %>%
        mutate(missing = sum(is.na(q1))) %>%
        filter(missing < 24) %>%
        ungroup() %>%
        select(-missing)
  

# Bootstrap ---------------------------------------------------------------
  # Bootstrap a sample of n = 100
  dl <- dl %>%
        spread(target, q1) %>%
        sample_n(size = 100, replace = TRUE) %>%
        mutate(participant = 1:n()) %>%
        gather("target", "q1", -participant)
  
  # Code target categories
  dl <- dl %>% mutate(
    category = case_when(
      target %in% 1:6   ~ 1L,
      target %in% 7:8   ~ 2L,
      target %in% 9:10  ~ 3L,
      target %in% 11:12 ~ 4L,
      target %in% 13:14 ~ 5L,
      target %in% 15:16 ~ 6L,
      target %in% 17:18 ~ 7L,
      target %in% 19:24 ~ 8L
    )
  )
  
  # Compose data list
  data_list <- list(
    N = 100*24,
    J = 100,
    K = 8,
    jj = dl$participant,
    kk = dl$category,
    y = dl$q1
  )

  
# Estimate ----------------------------------------------------------------
  m1_stan <- stan("models/q1_m1.stan", data = data_list, seed = seeds[1],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  

# Predict -----------------------------------------------------------------
  m1_stan %>%
    spread_draws(b_k[category]) %>% 
    summarise_at(vars(b_k), funs(median, sd))
  