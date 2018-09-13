rm(list = ls())

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes)


# Import ------------------------------------------------------------------
  # Observed proportions
  dl   <- read_rds("data/d1.rds") %>% 
          filter(!is.na(q1)) %>%
          mutate(ig_caste = recode_factor(
            ig_caste, "gm" = "GM", "obc" = "OBC", "scst" = "SC/ST", .ordered = TRUE
          )) %>%
          mutate(caste = recode_factor(
            caste, "GM" = "GM", "OBC" = "OBC", "SCST" = "SC/ST", "N/A" = NA_character_, .ordered = TRUE
          ))
  
  # Estimated probabilites
  post <- read_rds("models/q1_post.rds") %>%
          unnest(p_est) %>%
          mutate(ig_caste = recode_factor(
            ig_caste, "1" = "GM", "2" = "OBC", "3" = "SC/ST", .ordered = TRUE
          ))
  

# Summarise ---------------------------------------------------------------
  # Observed proportions
  d_obs <- dl %>% 
    group_by(target, category, lname, nation, religion, caste, ig_caste) %>%
    summarise(p_obs = mean(q1))
  
  # Estimated probabilities
  d_est <- post %>%
    rename(p_est = p_k) %>%
    group_by(model, category, ig_caste) %>%
    median_hdi(.width = c(.67, .89, .97))


# Plot --------------------------------------------------------------------
  p3 <- d_est %>%
    left_join(d_obs, ., by = c("category", "ig_caste")) %>%
    filter(model == 2) %>%
  ggplot(., aes(x = target, y = p_est, ymin = .lower, ymax = .upper)) +
    geom_ribbon(aes(group = .width), fill = "white", colour = NA) +
    geom_ribbon(aes(fill = ig_caste, group = 1-.width), colour = NA, alpha = 0.20) +
    geom_path(aes(colour = ig_caste)) +
    geom_point(aes(y = p_obs), shape = 3, colour = "grey60") +
    geom_point(aes(colour = ig_caste)) +
    geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5, 20.5), colour = "white", size = 2) +
    geom_vline(xintercept = 0.5, colour = "black") +
    scale_x_reverse(breaks = 1:24, minor_breaks = NULL) +
    scale_y_continuous(breaks = seq(0.00, 1.00, 0.50), minor_breaks = NULL) +
    scale_colour_viridis_d(begin = 0, end = 0.7) +
    scale_fill_viridis_d(begin = 0, end = 0.7) +
    coord_flip(xlim = c(0.5, 24.5), ylim = c(-0.1, 1.1), expand = FALSE) +
    facet_grid(. ~ ig_caste) +
    theme_grey(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(colour = "black", size = 10),
      axis.text.x = element_text(colour = "black", size = 10),
      axis.line.x = element_line(colour = "black"),
      axis.title = element_blank(),
      axis.ticks.x = element_line(colour = "black"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  

# Table -------------------------------------------------------------------
  t3 <- dl %>%
    distinct(target, nation, religion, caste) %>%
    mutate(
      `#` = ifelse(target < 10, paste(" ", target, sep = ""), target),
      nation = case_when(
        target %in% c(1, 5, 9, 13) ~ "Indian",
        target == 17 ~ "Nepali",
        target == 19 ~ "Sri Lankan",
        target == 21 ~ "Bangladeshi",
        TRUE ~ ""
      ),
      religion = case_when(
        target %in% seq(1, 21, 4) ~ religion,
        target == 19 ~ religion,
        TRUE ~ ""
      ),
      caste = case_when(
        target %in% seq(1, 13, 4) ~ as.character(caste),
        TRUE ~ ""
      )
    ) %>% 
    gather("column", "text", -target) %>% 
    mutate(
      column = recode_factor(
        column, "nation" = "Nation", "religion" = "Religion", "caste" = "Caste", "#" = "#", .ordered = TRUE
      ),
      hjust = ifelse(column == "#", 1, 0),
      position = ifelse(column == "#", 1.18, 0),
      .lower = 0,
      .upper = ifelse(column == "#", 1.18, 6)
    ) %>%
  ggplot(., aes(x = target, y = position, label = text)) +
    geom_linerange(aes(ymin = .lower, ymax = .upper), alpha = 0) +
    geom_text(aes(hjust = hjust), size = 10/.pt) +
    geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5, 20.5),
               colour = "black", size = 0.125) +
    geom_vline(xintercept = 0.5, colour = "black") +
    scale_x_reverse(expand = c(0, 0)) +
    coord_flip(xlim = c(0.5, 24.5)) +
    facet_grid(. ~ column, scales = "free_x", space = "free_x") +
    theme_grey(base_size = 10) +
    theme(
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(colour = "black", size = 10, hjust = 0),
      panel.spacing = unit(0.0, "lines"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.ticks = element_blank(),
      axis.text.x = element_text(colour = NA, size = 10),
      axis.text.y = element_blank()
    )


# Combine -----------------------------------------------------------------
  f3 <- gridExtra::grid.arrange(t3, p3, ncol = 2)
  
  # Export
  ggsave(f2, file = "figures/figure-3.pdf", width = 390/.pt, height = 160, units = "mm")
  