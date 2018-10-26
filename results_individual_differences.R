rm(list = ls())

# Notes -------------------------------------------------------------------
  # In total, I considered 5 * 4 = 20 items measuring intergroup contact. Across
  # items, less than 1.6% (n = 66) of individual responses were missing - though
  # about 13% of all participants failed to respond to at least one item. As 
  # item-by-item missingness was negligible, I used single imputation rather 
  # than (computationally demanding) full imputation. Before running Models 4 to 
  # 6, I imputed missing values from a multivariate normal distribution (see 
  # "impute_data.R"). For all subsequent analyses, I substituted the most likely
  # value (from the posterior distribution of the imputation model) for each 
  # missing value (in the original dataset).


# Library -----------------------------------------------------------------
  library(tidyverse); library(rstan); library(tidybayes); library(loo)

  # Stan options
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  seeds = c(6678002, 8970252,	9556204, 8660009)
  

# Functions ---------------------------------------------------------------
  inv_logit <- function(x) exp(x) / ( 1 + exp(x) )
  

# Import ------------------------------------------------------------------
  source("data/prepare_data_ic.R")
  source("data/prepare_data_sdo.R")

  # Compose data list
  data_list <- with(dl, list(
    # Numbers
    N = length(participant),
    J = n_distinct(participant),
    K = n_distinct(category),
    # Indices
    jj = participant,
    kk = category,
    # Vectors
    x_scst = ifelse(ig_caste == "scst", 1L, 0L),
    x_obc = ifelse(ig_caste == "obc", 1L, 0L),
    x_og = og,
    x_cq = replace_na(cq, -99), # Replace NAs for foreign targets 
    x_pc = replace_na(pc, -99), # Replace NAs for foreign targets
    x_nc = replace_na(nc, -99), # Replace NAs for foreign targets
    x_of = replace_na(of, -99), # Replace NAs for foreign targets
    x_gs = gs,
    x_sdo_d = sdo_d,
    x_sdo_e = sdo_e,
    y = q1
  ))


# Models ------------------------------------------------------------------
  # Model: M2 + Participant: Intergroup contact (4)
  m4_stan <- stan("models/q1_m4.stan", data = data_list, seed = seeds[1],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Posterior probabilities
  m4_stan %>% 
    gather_draws(b_cq, b_pc, b_nc, b_of) %>%
    mutate(.value = exp(.value)) %>%
    median_hdi(.width = .97)
  
  # Model: M2 + Participant: Intergroup contact (2)
  m5_stan <- stan("models/q1_m5.stan", data = data_list, seed = seeds[2],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Model: M5 + Participant: Intergroup contact (2) * Target: Category (4)
  m6_stan <- stan("models/q1_m6.stan", data = data_list, seed = seeds[3],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Model: M2 + Participant: SDO (2)
  m7_stan <- stan("models/q1_m7.stan", data = data_list, seed = seeds[4],
                  iter = 2000, warmup = 1000,
                  pars = c("alpha"), include = FALSE)
  
  # Posterior probabilities
  m7_stan %>%
    gather_draws(b_sdo_d, b_sdo_e) %>%
    mutate(.value = exp(.value)) %>%
    median_hdi(.width = .97)


# Compare -----------------------------------------------------------------
  # Compare models
  mc_fit <- list(m4_stan, m5_stan, m6_stan, m7_stan)

  # Merge with m0_stan to m3_stan
  mc_fit <- append(read_rds("models/q1_mc_fit.rds"), mc_fit)

  # Export
  write_rds(mc_fit, "models/q1_mc_fit.rds")


# Predict -----------------------------------------------------------------
  # Extract posterior samples
  post <- crossing(
      category = 1:4,
      ig_caste = 1:3,
      item = c("nc", "of"),
      response = seq(1, 5, 0.5)
    ) %>% filter(
      category != ig_caste,
      (item == "nc" & (response %% 1) == 0) | (item == "of")
    ) %>% 
    full_join(
      m5_stan %>% 
        spread_draws(l_k[category, ig_caste], b_nc, b_of) %>% 
        filter(category <= 4, category != ig_caste),
      by = c("category", "ig_caste")
    ) %>%
    left_join(ds, by = "item") %>%
    mutate(
      nc = ifelse(item == "of", 0, (response - mean) / sd),
      of = ifelse(item == "nc", 0, (response - mean) / sd),
      p_est = inv_logit(l_k + b_nc * nc + b_of * of)
    ) %>%
    select(.draw, category, ig_caste, item, response, p_est) %>%
    ungroup()
  
  # Export
  write_rds(post, "models/q1_ic_post.rds")
