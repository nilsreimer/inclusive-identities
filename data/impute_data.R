rm(list = ls())

# Notes -------------------------------------------------------------------


# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(tidybayes)

  # Stan options
  options(mc.cores = parallel::detectCores())
  seed = 70475734
  

# Functions ---------------------------------------------------------------
  # Standardise variables
  x_scale   <- function(v, y = v) (v - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
  x_unscale <- function(v, y) v * sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)


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
  ) %>% mutate_at(
    vars(q4:q60_8), x_scale
  )
  
  # Compose data list
  data_list <- list(
    # Matrix
    X = d_imp %>% select(-participant) %>% data.matrix(),
    # Numbers
    N_row = nrow(X),
    N_col = ncol(X),
    N_obs = length(X[!is.na(X)]),
    N_mis = length(X[is.na(X)]),
    # Indices
    ii_obs = matrix(1:length(X), nrow = nrow(X), ncol = ncol(X))[!is.na(X)],
    ii_mis = matrix(1:length(X), nrow = nrow(X), ncol = ncol(X))[is.na(X)],
    # Vectors
    x_obs   = X[!is.na(X)],
    x_mu    = c(rep(0, ncol(X)-2), mean(X[,"ig_obc"]), mean(X[,"ig_scst"])),
    x_sigma = c(rep(1, ncol(X)-2), sd(X[,"ig_obc"]), sd(X[,"ig_scst"]))
  )
  

# Model -------------------------------------------------------------------
  # Initial values
  estimates <- function(X, N_mis, perturb = FALSE){
    if(perturb) X <- X + rnorm(length(X), 0, 1)
    x_mu    <- as.numeric(apply(X, 2, mean, na.rm = TRUE))
    x_sigma <- as.numeric(apply(X, 2,   sd, na.rm = TRUE))
    x_mis   <- rnorm(N_mis, 0, 1)
    Lcorr   <- t(chol(cor(X, use = "pairwise.complete")))
    return(list(x_mu = x_mu, x_sigma = x_sigma, x_mis = x_mis, Lcorr = Lcorr))
  }
  inits <- function(chain_id){
    values <- estimates(data_list$X, data_list$N_mis, perturb = chain_id > 1)
    return(values)
  }

  # Model
  fit <- stan("data/impute_data_model.stan", data = data_list, init = inits,
              iter = 2000, warmup = 1000)

  
# Merge -------------------------------------------------------------------
  # Substitute best estimate in standardised data
  x_imp <- fit %>% spread_samples(x_mis[i]) %>% mode_qi()
  d_imp <- d_imp %>% 
    gather("item", "response", -participant) %>%
    mutate(i = 1:n()) %>%
    mutate(response = ifelse(i %in% data_list$ii_mis, x_imp$x_mis, response)) %>%
    select(-i) %>%
    spread(item, response) %>%
    select(-ig_obc, -ig_scst)
  
  # Transform imputed data to original response scale
  d_imp <- d_imp %>% 
    gather("item", "response", -participant) %>%
    left_join(
      d2 %>% select(names(d_imp)) %>% 
        gather("item", "response", -participant) %>% 
        group_by(item) %>% 
        summarise_at(vars(response), funs(mean, sd), na.rm = TRUE),
      by = "item"
    ) %>%
    mutate(response = response * sd + mean) %>%
    select(-mean, -sd) %>%
    spread(item, response)
  
  # Merge data set
  d2_impute <- d2 %>%
    select(-one_of(names(d_imp[,-1]))) %>%
    left_join(d_imp, by = "participant") %>%
    select(names(d2))


# Export ------------------------------------------------------------------
  write_rds(d2_impute, "data/d2_impute.rds")

  