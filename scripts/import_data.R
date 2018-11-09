rm(list = ls())

# Notes -------------------------------------------------------------------
  # See materials/materials.pdf (Appendix D) for items and answer options.

  # To protect participants, we do not provide the raw data (data/sessions/)
  # with verbatim comments; this R script is thus not functional. Please import 
  # data/d1_raw.rds and data/d2_raw.rds directly.


# Library -----------------------------------------------------------------
  library(tidyverse)


# Functions ---------------------------------------------------------------
  clean <- function(v, type, range = FALSE) {
    if (is_character(v)) {
      if (type == "n") {
        v <- as.numeric(v)
      } else if (type == "i") {
        v <- as.integer(v)
      } else if (type == "c") {
        v <- as.character(v)
      }
      if (length(range) > 1) {
        v <- ifelse(v < min(range) | v > max(range), NA, v)
      }
    } else {
      warning("Input is not a character vector.")
    }
    return(v)
  }

  
# Import ------------------------------------------------------------------
  # Import and merge csv files
  for (session in 1:8) {
    if (session == 1) {
      d1 <- read_csv(
        paste("data/sessions/", LETTERS[session], "_DL_PARTI.csv", sep=""), 
        col_types=cols(.default=col_character())
      )  
      d2 <- read_csv(
        paste("data/sessions/", LETTERS[session], "_WL_PARTII.csv", sep=""), 
        col_types=cols(.default=col_character())
      ) 
    } else {
      d1 <- bind_rows(d1, read_csv(
        paste("data/sessions/", LETTERS[session], "_DL_PARTI.csv", sep=""), 
        col_types=cols(.default=col_character()))
      )
      d2 <- bind_rows(d2, read_csv(
        paste("data/sessions/", LETTERS[session], "_WL_PARTII.csv", sep=""), 
        col_types=cols(.default=col_character()))
      )
    }
  }
  
  # Lower-case variable names
  d1 <- d1 %>% rename_all(tolower)
  d2 <- d2 %>% rename_all(tolower)
  
  # Correct coding error in Session H
  d1 <- bind_rows(
    d1 %>% filter(session != "H"),
    d1 %>% filter(session == "H") %>% 
           mutate(participant = as.character(rep(1:54, each = 24)))
  )
  d2 <- bind_rows(
    d2 %>% filter(session != "H"),
    d2 %>% filter(session == "H") %>% mutate(participant = as.character(1:54))
  )


# Clean -------------------------------------------------------------------
  # Part 1 (TCCT)
  d1 <- d1 %>% mutate(
    session = clean(session, "c"),
    participant = clean(participant, "i"),
    target = clean(target, "i", 1:24),
    q1 = clean(q1, "i", 0:1),
    q2 = clean(q2, "i", 1:7),
    q3 = clean(q3, "i", 0:100)
  )
  
  # Part 2 (Questionnaire)
  d2 <- d2 %>% 
    mutate(
      session = clean(session, "c"), 
      participant = clean(participant, "i")
    ) %>% 
    mutate_at(vars(q4:q23), ~clean(., "i", 1:5)) %>%
    mutate_at(vars(q24_1:q24_7), ~clean(., "i", 1:7)) %>%
    mutate_at(vars(q25:q29), ~clean(., "i", 1:5)) %>%
    mutate_at(vars(q30_1:q31_5), ~clean(., "i", 1:5)) %>%
    mutate(
      q32 = clean(q32, "i", 1:7),
      q33 = clean(q33, "i", 1:3),
      q34 = clean(q34, "i", 1:10),
      q35_1 = clean(q35_1, "i", 0:1),
      q35_2 = clean(q35_2, "c"),
      q36 = clean(q36, "i", 1:7),
      q37 = clean(q37, "i", 1:5)
    ) %>%
    mutate_at(vars(q38_1:q38_3), ~clean(., "i", 1:7)) %>%
    mutate(
      q39 = clean(q39, "i", 1:5),
      q40 = clean(q40, "i", 1:5),
      q41 = clean(q41, "i", 0:100),
      q42 = clean(q42, "i", 0:100)
    ) %>%
    mutate_at(vars(q43:q52), ~clean(., "i", 1:7)) %>%
    mutate(
      q53 = clean(q53, "i", 1:5),
      q54 = clean(q54, "i", 1:5),
      q55 = clean(q55, "i", 0:100),
      q56 = clean(q56, "i", 0:100),
      q57 = clean(q57, "i", 0:100)
    ) %>%
    mutate_at(vars(q58_1:q60_8), ~clean(., "i", 1:7))
  

# Match -------------------------------------------------------------------
  # Import session-wise target information
  for (session in 1:8) {
    if (session == 1) {
      p1 <- read_csv(paste("materials/sessions/SESSION_", LETTERS[session], 
                           ".csv", sep = "")) %>% 
            rename_all(tolower) %>%
            mutate(session = LETTERS[session])
    } else {
      p1 <- read_csv(paste("materials/sessions/SESSION_", LETTERS[session], 
                           ".csv", sep = "")) %>% 
            rename_all(tolower) %>%
            mutate(session = LETTERS[session]) %>%
            bind_rows(p1, .)
    }
  }
  
  # Match with Part 1 data
  d1 <- left_join(d1, p1, by = c("session", "target")) %>%
    select(session, participant, target, index:caste, q1:q3, comments)
  
  
# Transform ---------------------------------------------------------------
  # Assign each participant an unambigious id
  d1 <- d1 %>% mutate(participant = paste(session, participant, sep = ""))
  d2 <- d2 %>% mutate(participant = paste(session, participant, sep = ""))
  
  # Rename target and index to make column names clearer.
  d1 <- d1 %>% rename(order = target, target = index)


# Export ------------------------------------------------------------------
  write_rds(d1 %>% select(-comments), "data/d1_raw.rds")
  write_rds(d2 %>% select(-feedback, -comments), "data/d2_raw.rds")
