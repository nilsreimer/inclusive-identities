rm(list = ls())

# Notes -------------------------------------------------------------------
  # This script imputes missing responses (n = 66) to the intergroup contact
  # variables, and returns a data frame with both observed responses and 
  # draws from the posterior distribution for missing responses. 

# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(tidybayes)

  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seed = 70475734


# Import ------------------------------------------------------------------
  d2 <- read_rds("data/d2.rds") %>%
        select(participant, ig_caste, q4:q23)
  
  # Prepare data (long format)
  dl <- d2 %>%
    gather("item", "response", q4:q23) %>%
    mutate(item = str_extract(item, "[0-9]+") %>% as.integer()) %>% 
    arrange(participant, item) %>%
    mutate(ii = 1:n()) %>%
    select(ii, participant, item, response)
  
  # Prepare data (wide format)
  dw <- d2 %>%
    mutate(
      x_obc = ifelse(ig_caste == "obc", 1L, 0L),
      x_scst = ifelse(ig_caste == "scst", 1L, 0L)
    ) %>%
    select(participant, x_obc, x_scst) %>%
    arrange(participant)
  
  # Summarise
  ds <- dl %>% 
    group_by(item) %>% 
    summarise_at(vars(response), funs(mean, sd), na.rm = TRUE)
  
  # Compose data list
  data_list <- list(
    # Numbers
    N_row = n_distinct(dl$participant),
    N_col = n_distinct(dl$item),
    N_mis = sum(is.na(dl$response)),
    N_obs = sum(!is.na(dl$response)),
    # Indices
    ii_mis  = dl$ii[is.na(dl$response)],
    ii_obs  = dl$ii[!is.na(dl$response)],
    ii_item = dl$item[!is.na(dl$response)] - 3,
    # Vectors
    x_obc = dw$x_obc,
    x_scst = dw$x_scst,
    x_mean = ds$mean,
    x_sd = ds$sd,
    x_obs = dl$response[!is.na(dl$response)]
  )
  
# Model -------------------------------------------------------------------
  # Initial values
  X <- dl %>% 
       select(-ii) %>% 
       spread(item, response) %>% 
       select(-participant) %>% 
       data.matrix()
  estimates <- function(X, data_list, perturb = FALSE){
    if(perturb) X <- X + rnorm(length(X), 0, 1)
    list(
      b_0     = rnorm(data_list$N_col, 0, 1),
      b_obc   = rnorm(data_list$N_col, 0, 1),
      b_scst  = rnorm(data_list$N_col, 0, 1),
      x_sigma = rep(1, data_list$N_col),
      x_mis   = rnorm(data_list$N_mis, 0, 1),
      Lcorr   = t(chol(cor(X, use = "pairwise.complete")))
    ) %>% return()
  }
  inits <- function(chain_id){
    values <- estimates(X, data_list, perturb = chain_id > 1)
    return(values)
  }
  
  # Model
  fit <- stan("models/impute_data.stan", data = data_list, seed = seed,
              iter = 2000, warmup = 1000,
              pars = c("X", "Mu", "Lcorr"), include = FALSE)


# Prepare -----------------------------------------------------------------
  # Extract posterior draws
  d_imp <- fit %>% 
    spread_draws(x_mis[ii]) %>%
    ungroup() %>% 
    mutate(ii = map_int(ii, ~data_list$ii_mis[.])) %>%
    rename(x_imp = x_mis) %>%
    select(.draw, ii, x_imp)
  
  # Rescale posterior draws
  d_imp <- d_imp %>%
    left_join(dl %>% select(ii, item), by = "ii") %>%
    left_join(ds, by = "item") %>%
    mutate(x_imp = x_imp * sd + mean) %>%
    select(.draw, ii, x_imp)
  
  # Code contact variables
  d_obs <- dl %>%
    mutate(
      category = case_when(
        item <=  8 ~ 3L,
        item <= 13 ~ 2L,
        item <= 18 ~ 1L,
        item <= 23 ~ 4L
      ),
      item = case_when(
        item %in% seq(4, 19, 5) ~ "cq",
        item %in% seq(5, 20, 5) ~ "pc",
        item %in% seq(6, 21, 5) ~ "nc",
        TRUE ~ "of"
      ),
      x_obs = response,
      missing = is.na(response)
    ) %>%
    select(ii, participant, category, item, x_obs, missing)
  
  # Merge with observed responses
  d_mis <- d_obs %>%
    group_by(participant, category, item) %>%
    filter(mean(missing) > 0) %>%
    expand(.draw = 1:4000, ii) %>%
    ungroup() %>%
    left_join(d_imp, by = c(".draw", "ii")) %>%
    left_join(d_obs, by = c("participant", "category", "item", "ii")) %>%
    mutate(response = ifelse(missing, x_imp, x_obs)) %>%
    group_by(participant, category, item, .draw) %>%
    summarise(response = mean(response)) %>%
    ungroup() %>%
    select(-.draw) %>%
    mutate(missing = TRUE)
  
  # Observed responses
  d_obs <- d_obs %>%
    group_by(participant, category, item) %>%
    filter(mean(missing) == 0) %>%
    summarise(response = mean(x_obs)) %>%
    ungroup() %>%
    mutate(missing = FALSE)
  
  # Merge imputed and observed responses
  d2_impute <- bind_rows(d_obs, d_mis) %>%
    left_join(
      d_obs %>% group_by(item) %>% summarise_at(vars(response), funs(mean, sd)),
      by = "item"
    ) %>%
    mutate(z_response = (response - mean) / sd) %>%
    select(participant, category, item, response, missing, z_response, mean, sd)
  

# Export ------------------------------------------------------------------
  d2_impute %>% 
    mutate(item = ordered(item, levels = c("cq", "pc", "nc", "of"))) %>%
    arrange(participant, category, item) %>%
    write_rds("data/d2_impute.rds")
