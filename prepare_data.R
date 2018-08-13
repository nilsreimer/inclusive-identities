rm(list = ls())

# Notes -------------------------------------------------------------------
  #  See materials/script.pdf (Appendix D) for items and answer options.


# Library -----------------------------------------------------------------
  library(tidyverse)


# Functions ---------------------------------------------------------------


  
# Import ------------------------------------------------------------------
  d1 <- read_rds("data/d1_raw.rds")
  d2 <- read_rds("data/d2_raw.rds")

  
# Participants ------------------------------------------------------------
  # Count religion, caste (after exclusion)
  d2 %>% count(q36)
  d2 %>% count(q37)
  
  # Exclude participants
  d2 <- d2 %>% filter(q36 != 2, q37 != 5)
  d1 <- d1 %>% filter(participant %in% d2$participant)
  
  # Count religion, caste (after exclusion)
  d2 %>% count(q36)
  d2 %>% count(q37)
  

# Transform ---------------------------------------------------------------
  # Code ingroup caste
  d2 <- d2 %>% mutate(
    ig_caste = recode(q37, "1" = "scst", "2" = "scst", "3" = "obc", "4" = "gm")
  ) %>% select(
    session, participant, ig_caste, everything()
  )
  
  # Code target categories
  d1 <- d1 %>% mutate(
    category = case_when(
      target <= 4 ~ 1L,
      target <= 8 ~ 2L,
      target <=12 ~ 3L,
      target <=16 ~ 4L,
      target <=20 ~ 5L,
      target <=24 ~ 6L
    )
  ) %>% left_join(d2 %>% select(participant, ig_caste)) %>% select(
    session:target, category, lname:caste, ig_caste, everything()
  )

  # Convert participant id to integer
  d1 <- d1 %>% mutate(participant = as.integer(ordered(participant)))
  d2 <- d2 %>% mutate(participant = as.integer(ordered(participant)))
  

# Export ------------------------------------------------------------------
  write_rds(d1, "data/d1.rds")
  write_rds(d2, "data/d2.rds")
  
  
# Prepare -----------------------------------------------------------------
  # Impute missing data (predictor variables)
  # source("data/impute_data.R") 
  # This takes a long time; import read_rds("data/d2_impute.rds") instead.
  
  # Prepare data for analyses
  source("data/prepare_data_ic.R")
  source("data/prepare_data_sdo.R")
  source("data/prepare_data_it.R")
  
  # Compile data for analyses
  # d1
  # d1_pred
  # d2_it
  # d2_ps

# Export ------------------------------------------------------------------
  write_rds(d1, "data/d1_with_predictors.rds")
  write_rds(d2_it, "data/d2_intergroup_threat.rds")
  
  # Export to other formats
  # write_csv(d1, "data/export/d1.csv")
  # write_csv(d1_pred, "data/export/d1_with_predictors.csv")
  # write_csv(d2, "data/export/d2.csv")
