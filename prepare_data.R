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
  # Count religion, caste (before exclusion)
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
    ig_caste = recode(q37, "1" = "SCST", "2" = "SCST", "3" = "OBC", "4" = "GM")
  ) %>% select(
    session, participant, ig_caste, everything()
  )
  
  # Code target categories
  d1 <- d1 %>% mutate(
    category = case_when(
      index <= 4 ~ 1L,
      index <= 8 ~ 2L,
      index <=12 ~ 3L,
      index <=16 ~ 4L,
      index <=20 ~ 5L,
      index <=24 ~ 6L
    )
  ) %>% left_join(d2 %>% select(participant, ig_caste)) %>% select(
    session:index, category, lname:caste, ig_caste, everything()
  )


# Export ------------------------------------------------------------------
  write_rds(d1, "data/d1.rds")
  write_rds(d2, "data/d2.rds")
  

# Prepare -----------------------------------------------------------------
  # Individual differences
  source("data/prepare_data_ic.R")
  source("data/prepare_data_sdo.R")
  
  # Consequences
  source("data/prepare_data_it.R")
