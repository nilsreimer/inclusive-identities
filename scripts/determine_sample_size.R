rm(list = ls())

# Notes -------------------------------------------------------------------

  #########################################################################
  # This R script is included for reference only as we are not able to    #
  # make the underlying data available. For details, see Appendix B.      #
  #########################################################################

# Library -----------------------------------------------------------------
  
  # Load packages
  library(tidyverse); library(rstan); library(tidybayes)
  
  # Stan options
  n_cores <- 8L
  options(mc.cores = n_cores)
  rstan_options(auto_write = TRUE)
  seeds = c(35702546)
  
  # Link function
  inv_logit <- function(x) exp(x) / (1 + exp(x))


# Prepare -----------------------------------------------------------------

  # Import data in wide format
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
    select(-missing) %>% 
    mutate(participant = as.integer(factor(participant)))
  
  # Code target categories
  dl <- dl %>% 
    mutate(
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
    ) %>% 
    select(participant, target, category, q1)


# Simulate ----------------------------------------------------------------

  # Compose data list
  data_list <- with(dl, list(
    N = length(participant),
    J = n_distinct(participant),
    K = n_distinct(category),
    jj = participant,
    kk = category,
    y = q1
  ))
  
  # Set seed for simulation
  set.seed(35702546)
  
  # Simulate data
  dl_sim <- dl %>% 
    distinct(category, target) %>% 
    crossing(
      nesting(
        participant = 1:300,
        x = rnorm(300, 0, 1),
        group = rep(1:3, each = 100)
      )
    ) %>% 
    arrange(participant, target) %>% 
    rowid_to_column("ii")
  
  # Append simulation parameters
  data_list <- with(dl_sim, list(
    N_sim  = length(participant),
    J_sim  = n_distinct(participant),
    K_sim  = n_distinct(category),
    jj_sim = participant,
    kk_sim = category,
    x_sim = x,
    x_g1 = if_else(group == 1L, 1L, 0L),
    x_g2 = if_else(group == 2L, 1L, 0L)
  )) %>% append(data_list, .)
  
  # Run model
  sim <- stan(
    "models/simulate_samples.stan", 
    data = data_list, 
    iter = 500 + 4000/n_cores, warmup = 500, 
    chains = n_cores,
    control = list(adapt_delta = 0.99),
    seed = seeds[1],
    pars = c("alpha"), include = FALSE
  )
  

# Estimate ----------------------------------------------------------------

  # Combine simulations
  dl_est <- dl_sim %>% 
    crossing(.draw = 1:4000) %>% 
    arrange(.draw, ii) %>%
    left_join(
      as.data.frame(sim, "y_sim") %>% 
        rowid_to_column(".draw") %>% 
        pivot_longer(
          -.draw,
          names_to = "ii", 
          names_pattern = "y_sim\\[([0-9]*)\\]",
          names_ptypes = list(ii = integer()),
          values_to = "y_sim",
          values_ptypes = list(y_sim = integer())
        ),
      by = c(".draw", "ii")
    )
    
  # Compose data list
  make_dlist <- function(d) {
    dl_est %>% 
      filter(.draw == d) %>% 
      with(., list(
        N = length(participant),
        J = n_distinct(participant),
        K = n_distinct(category),
        jj = participant,
        kk = category,
        y = y_sim,
        x_sim = x,
        x_g1 = if_else(group == 1L, 1L, 0L),
        x_g2 = if_else(group == 2L, 1L, 0L)
      ))
  }
  
  # Select samples
  samples <- sample.int(4000, 10)
  
  # Run models
  start_time <- Sys.time()
  for (i in 1:length(samples)) {
    print(glue::glue("[{i}/{length(samples)}], {lubridate::as.duration(Sys.time() - start_time)}"))
    fit <- stan(
      "models/estimate_samples.stan", 
      data = make_dlist(samples[i]), 
      iter = 500 + 2000/n_cores, warmup = 500, 
      chains = n_cores,
      pars = c("alpha", "b_j"), include = FALSE
    )
    if (i == 1) {
      results <- list(fit)
    } else {
      results <- append(results, fit)
    }
  }
  
  # Compile results
  sim_fit <- tibble(sample = 1:length(samples), fit = results)

  # Save results
  write_rds(sim_fit, "results/sim_fit.rds")
  

# Analyze -----------------------------------------------------------------

  # Load results
  sim_fit <- read_rds("results/sim_fit.rds")
  
  # Calculate and summarize SDs of posterior probabilities for predictors
  sim_fit %>% 
    mutate(
      b_x = map_dbl(
        fit,
        ~spread_draws(.x, b_x) %>% pull(b_x) %>% sd()
      ),
      b_g = map(
        fit,
        ~ spread_draws(
          .x,
          b_g1, b_g1_k[kk], sigma_g1_k, 
          b_g2, b_g2_k[kk], sigma_g2_k
        ) %>% 
          mutate(
            b_g1 = b_g1 + b_g1_k*sigma_g1_k,
            b_g2 = b_g2 + b_g2_k*sigma_g2_k
          ) %>% 
          select(.draw, kk, b_g1, b_g2) %>% 
          ungroup() %>% 
          pivot_wider(
            names_from = kk,
            values_from = c(b_g1, b_g2)
          ) %>% 
          summarise_at(vars(-.draw), sd)
      )
    ) %>% 
    unnest(b_g) %>% 
    pivot_longer(
      -sample:-fit,
      names_to = "parameter",
      values_to = "sd"
    ) %>% 
    group_by(parameter) %>% 
    summarise_at(
      vars(sd),
      list("min" = min, "max" = max, "sd" = mean)
    ) %>% 
    mutate(
      p_sd = inv_logit(0.96 + sd) - inv_logit(0.96)
    ) %>% 
    # group_by(str_detect(parameter, "b_g")) %>% summarise_at(vars(sd, min, max, p_sd), mean) %>% 
    mutate_at(vars(sd, min, max), round, digits = 2) %>% 
    mutate_at(vars(p_sd), round, digits = 3)
    
  