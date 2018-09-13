# Import Part 1 data (long)
  dl <- read_rds("data/d1.rds") %>%
        filter(!is.na(q1)) %>%
        mutate(og = case_when(
          category == 1 & ig_caste == "gm" ~ 0L,
          category == 2 & ig_caste == "obc" ~ 0L,
          category == 3 & ig_caste == "scst" ~ 0L,
          TRUE ~ 1L
        )) %>%
        arrange(participant, target)

# Import Part 2 data (intergroup contact)
  d2_imp <- read_rds("data/d2_impute.rds") %>%
            group_by(participant, category, item) %>%
            summarise(x_imp = mean(z_response)) %>%
            spread(item, x_imp)

# Merge Part 1 and Part 2 data (intergroup contact)
  dl <- left_join(dl, d2_imp, by = c("participant", "category"))
  
# Summaries used to standardise contact variables
  ds <- read_rds("data/d2_impute.rds") %>% distinct(item, mean, sd)
  
# Export
  write_rds(dl, "dl_ic.rds")
  rm("d2_imp")

