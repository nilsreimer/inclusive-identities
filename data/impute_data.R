rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(tidybayes)

  # Stan options
  options(mc.cores = parallel::detectCores())
  seed = 70475734


# Import ------------------------------------------------------------------
  d2 <- read_rds("data/d2.rds")

  
# Prepare -----------------------------------------------------------------
  # Prepare variables for imputation
  d_imp <- d2 %>% mutate(
    ig_gm   = ifelse(ig_caste ==   "gm", 1L, 0L),
    ig_obc  = ifelse(ig_caste ==  "obc", 1L, 0L),
    ig_scst = ifelse(ig_caste == "scst", 1L, 0L)
  ) %>% select(
    participant, ig_obc, ig_scst, q4:q23, q60_1:q60_8
  ) %>% gather("item", "response", -participant)
  
  # Standardise variables
  d_imp <- d_imp %>% left_join(
    d_imp %>% 
      group_by(item) %>% 
      summarise_at(vars(response), funs(mean, sd), na.rm = TRUE), 
    by = "item"
  ) %>% mutate(
    i = 1:n(),
    x = (response - mean) / sd
  )
  
  # Compose data list
  data_list <- list(
    # Numbers
    N_row = n_distinct(d_imp$participant),
    N_col = n_distinct(d_imp$item),
    N_mis = sum(is.na(d_imp$x)),
    N_obs = sum(!is.na(d_imp$x)),
    # Indices
    ii_mis = d_imp$i[is.na(d_imp$x)],
    ii_obs = d_imp$i[!is.na(d_imp$x)],
    # Vectors
    x_obs   = d_imp$x[!is.na(d_imp$x)]
  )
  

# Model -------------------------------------------------------------------
  # Initial values
  X <- d_imp %>% select(participant, item, x) %>% spread(item, x) %>% select(-participant) %>% data.matrix()
  estimates <- function(X, N_mis, perturb = FALSE){
    if(perturb) X <- X + rnorm(length(X), 0, 1)
    x_mu    <- rep(0, 30)
    x_sigma <- rep(1, 30)
    x_mis   <- rnorm(N_mis, 0, 1)
    Lcorr   <- t(chol(cor(X, use = "pairwise.complete")))
    return(list(x_mu = x_mu, x_sigma = x_sigma, x_mis = x_mis, Lcorr = Lcorr))
  }
  inits <- function(chain_id){
    values <- estimates(X, data_list$N_mis, perturb = chain_id > 1)
    return(values)
  }

  # Model
  fit <- stan("data/impute_data_model.stan", data = data_list, init = inits,
              iter = 2000, warmup = 1000, seed = seed)

  
# Merge -------------------------------------------------------------------
  # Substitute best estimate in standardised data
  x_imp <- fit %>% spread_samples(x_mis[i]) %>% mode_qi() %>% ungroup()
  d_imp <- d_imp %>%
    left_join(
      x_imp %>% mutate(i = data_list$ii_mis), 
      by = "i"
    ) %>%
    mutate(
      x = ifelse(is.na(x), x_mis, x),
      x = x * sd + mean,
      response = ifelse(is.na(response), x, response)
    ) %>% select(participant, item, response) %>% 
    spread(item, response) %>%
    select(-starts_with("ig")) %>% 
    select(participant, q4, q5, q6, q7, q8, q9, everything())

# Export ------------------------------------------------------------------
  write_rds(d_imp, "data/d2_impute.rds")
  