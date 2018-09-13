rm(list = ls())

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes)


# Import ------------------------------------------------------------------
  post <- read_rds("models/q2q3_post.rds") %>%
    gather("standardized", "m_est", m_est:m_est_z) %>%
    mutate(standardized = (standardized == "m_est_z")) %>%
    mutate(
      ig_caste = recode_factor(
        ig_caste, 
        "1" = "GM",
        "2" = "OBC",
        "3" = "SC/ST",
        .ordered = TRUE
      ),
      category = recode_factor(
        category,
        "6" = "Foreign, Muslim",
        "5" = "Foreign, Hindu",
        "4" = "Indian, Muslim, OBC",
        "3" = "Indian, Hindu, SC/ST",
        "2" = "Indian, Hindu, OBC",
        "1" = "Indian, Hindu, GM",
        .ordered = TRUE
      ),
      x_q1 = recode_factor(
        x_q1,
        "1" = "us",
        "0" = "not us",
        .ordered = TRUE
      )
    ) %>%
    group_by(category, ig_caste, x_q1, x_obc, x_scst, item, standardized)
  

# Figure 5 ----------------------------------------------------------------
  post %>%
    filter(!standardized, item == "q2") %>%
  ggplot(., aes(x = m_est, y = category, group = x_q1)) +
    geom_halfeyeh(
      aes(shape = x_q1),
      colour = NA, fill = "white",
      point_interval = median_hdi, .width = 0.97,
      size = 0
    ) +
    geom_halfeyeh(
      aes(shape = x_q1, colour = ig_caste, fill = ig_caste), 
      point_interval = median_hdi, .width = 0.97,
      size = 1
    ) +
    scale_x_continuous(breaks = seq(1, 7, 3), minor_breaks = seq(1, 7, 1)) +
    scale_fill_viridis_d(begin = 0, end = 0.7, alpha = 2/5) +
    scale_colour_viridis_d(begin = 0, end = 0.7) +
    coord_cartesian(xlim = c(0.5, 7.5), ylim = c(0.9, 7), expand = FALSE) +
    facet_grid(. ~ ig_caste) + labs(
      x = "Social distance"
    ) +
    theme_grey(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.y = element_blank(),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black", hjust = 0, vjust = 0),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_line(size = 1)
    )
  
  # Export
  ggsave(file = "figures/figure-5.pdf", width = 390/.pt, height = 80, units = "mm")
  

# Figure 6 ----------------------------------------------------------------
  post %>%
    filter(!standardized, item == "q3") %>%
  ggplot(., aes(x = m_est, y = category, group = x_q1)) +
    geom_halfeyeh(
      aes(shape = x_q1),
      colour = NA, fill = "white",
      point_interval = median_hdi, .width = 0.97,
      size = 0
    ) +
    geom_halfeyeh(
      aes(shape = x_q1, colour = ig_caste, fill = ig_caste), 
      point_interval = median_hdi, .width = 0.97,
      size = 1
    ) +
    scale_x_continuous(breaks = seq(0, 100, 50), minor_breaks = seq(0, 100, 25)) +
    scale_fill_viridis_d(begin = 0, end = 0.7, alpha = 2/5) +
    scale_colour_viridis_d(begin = 0, end = 0.7) +
    coord_cartesian(xlim = c(-5, 106), ylim = c(0.9, 7), expand = FALSE) +
    facet_grid(. ~ ig_caste) + labs(
      x = "Feeling thermometer"
    ) +
    theme_grey(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.y = element_blank(),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black", hjust = 0, vjust = 0),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_line(size = 1)
    )
  
  # Export
  ggsave(file = "figures/figure-6.pdf", width = 390/.pt, height = 80, units = "mm")
  
  