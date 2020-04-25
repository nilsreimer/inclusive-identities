# Load packages -----------------------------------------------------------
  library(lavaan)


# Import ------------------------------------------------------------------
  d2 <- read_rds("data/d2.rds") %>%
        select(participant, starts_with("q60")) %>%
        arrange(participant) %>%
        rename_at(vars(starts_with("q60")), ~str_replace(., "q60", "sdo")) %>%
        mutate_at(vars(sdo_3, sdo_4, sdo_7, sdo_8), ~ 8 - .)


# Models ------------------------------------------------------------------
  m4_model <- '
  # Loadings
  sdo_d =~ sdo_1 + sdo_2 + sdo_3 + sdo_4
  sdo_e =~ sdo_5 + sdo_6 + sdo_7 + sdo_8
  pro_t =~ sdo_1 + sdo_2 + sdo_5 + sdo_6
  con_t =~ sdo_3 + sdo_4 + sdo_7 + sdo_8

  # Covariances
  sdo_d ~~ 0*pro_t
  sdo_d ~~ 0*con_t
  sdo_e ~~ 0*pro_t
  sdo_e ~~ 0*con_t

  # Variances
  sdo_d ~~ 1*sdo_d
  sdo_e ~~ 1*sdo_e
  pro_t ~~ 1*pro_t
  con_t ~~ 1*con_t
  '
  m4_fit <- cfa(m4_model, estimator = "mlr", missing = "fiml", data = d2)


# Export ------------------------------------------------------------------
  dl <- lavPredict(m4_fit) %>% 
    matrix(., nrow = 302, ncol = 4, dimnames = list(NULL, c("sdo_d", "sdo_e", "pro_t", "con_t"))) %>%
    as_tibble() %>% 
    bind_cols(d2, .) %>%
    select(participant, sdo_d, sdo_e) %>% 
    left_join(dl, ., "participant")
  
  # Transform
  dl <- dl %>% mutate(
    gs = case_when(
      ig_caste == "gm"   & category %in% 2:6 ~ 1L,
      ig_caste == "obc"  & category %in% 3:6 ~ 1L,
      ig_caste == "scst" & category %in% 4:6 ~ 1L,
      TRUE ~ 0L
    )
  ) %>% select(
    session:og, gs, sdo_d:sdo_e
  )


# Remove ------------------------------------------------------------------
  detach(package:lavaan)
  rm("d2", "m4_fit", "m4_model")
