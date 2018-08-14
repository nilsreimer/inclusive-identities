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


# Compare -----------------------------------------------------------------
  # Compare models
  mc_fit <- list(m0_stan, m1_stan, m2_stan, m3_stan)
  mc_llk <- map(mc_fit, extract_log_lik, merge_chains = FALSE)
  mc_loo <- map(mc_llk, ~loo(., r_eff = relative_eff(exp(.)), cores = 4))

  compare(x = mc_loo)
  
  crossing(model1 = 0:3, model2 = 0:3) %>% mutate(
    elpd_diff = map2_dbl(model1, model2, ~compare(mc_loo[[.y+1]], mc_loo[[.x+1]])["elpd_diff"]),
    se = map2_dbl(model1, model2, ~compare(mc_loo[[.y+1]], mc_loo[[.x+1]])["se"]),
    ratio = elpd_diff / se
  ) %>% ggplot(., aes(x = model2, y = model1)) +
    geom_tile(aes(fill = ratio)) +
    geom_text(aes(label = round(ratio, 1))) +
    scale_y_reverse() +
    scale_fill_gradient2(na.value = NA) +
    coord_fixed() +
    theme_minimal()
  
  tibble(
    model1 = 0:3,
    model2 = c(0L, 0L, 1L, 2L)
  ) %>% mutate(
    elpd = map_dbl(model1, ~mc_loo[[.+1]]$estimates["elpd_loo", "Estimate"]),
    se_elpd = map_dbl(model1, ~mc_loo[[.+1]]$estimates["elpd_loo", "SE"]),
    d_elpd = map2_dbl(model1, model2, ~compare(mc_loo[[.y+1]], mc_loo[[.x+1]])["elpd_diff"]),
    se_d_elpd = map2_dbl(model1, model2, ~compare(mc_loo[[.y+1]], mc_loo[[.x+1]])["se"]),
    r_d_elpd = d_elpd / se_d_elpd
  ) %>% mutate_if(is.double, round, digits = 1)
  
  # Export
  write_rds(mc_fit, "models/q1_mc_fit.rds")
  write_rds(mc_loo, "models/q1_mc_loo.rds")


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
  
  # post %>%
  #   unnest(p_est) %>%
  #   mutate(model = ordered(model)) %>%
  #   group_by(model, category, ig_caste) %>%
  #   # median_hdi(.width = .80) %>%
  # ggplot(., aes(y = category, x = p_k)) +
  #   geom_halfeyeh(aes(group = model, fill = model, colour = model)) +
  #   scale_fill_viridis_d(alpha = 0.25) +
  #   scale_y_reverse(breaks = 1:6) +
  #   facet_grid(. ~ ig_caste)
