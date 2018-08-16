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
  group_by(participant, category, item, missing) %>%
  summarise(x_imp = mean(z_response), x_imp_sd = sd(z_response)) %>%
  mutate(x_imp_sd = ifelse(missing, x_imp_sd, 0)) %>%
  ungroup()

# Merge Part 1 and Part 2 data (intergroup contact)
dl <- d2_imp %>%
  select(participant, category, item, x_imp) %>%
  spread(item, x_imp) %>%
  left_join(dl, ., by = c("participant", "category"))

# Prepare imputed data (intergroup contact)
dl_imp <- dl %>%
  select(participant, target, category, cq:of) %>%
  gather(item, x_imp, cq:of) %>%
  mutate(
    ii = 1:n(),
    item = ordered(item, levels = c("cq", "pc", "nc", "of"))
  ) %>% 
  left_join(
    d2_imp %>% select(-x_imp), 
    by = c("participant", "category", "item")
  ) %>% 
  select(ii, participant:item, missing, x_imp, x_imp_sd) %>%
  replace_na(list(x_imp = -99, missing = FALSE, x_imp_sd = 0))