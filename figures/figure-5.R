rm(list = ls())

# Library -----------------------------------------------------------------

  # Load packages
  library(tidyverse); library(patchwork); library(tidybayes)


# Prepare -----------------------------------------------------------------

  # Import posterior draws
  post <- read_rds("results/q2q3_q1_post.rds")

  # Recode variables
  post <- post %>% 
    mutate(
      category = recode_factor(
        category,
        "1" = "Indian, Hindu, GM",
        "2" = "Indian, Hindu, OBC",
        "3" = "Indian, Hindu, SC/ST",
        "4" = "Indian, Muslim, OBC",
        "5" = "Foreign, Hindu",
        "6" = "Foreign, Muslim"
      ),
      ig_caste = factor(ig_caste, levels = c("GM", "OBC", "SC/ST"))
    )


# Visualize ---------------------------------------------------------------

  # Figure 5a
  f5a <- post %>% 
    filter(outcome == "q2") %>% 
  ggplot(., aes(x = m_est, y = fct_rev(category), group = q1)) +
    stat_halfeyeh(
      aes(shape = factor(q1)),
      .width = .97, size = 0.25, fatten_point = 1,
      colour = "white", fill = "white"
    ) +
    stat_halfeyeh(
      aes(colour = ig_caste, fill = ig_caste, shape = factor(q1)),
      .width = .97, size = 0.25, fatten_point = 1,
      slab_alpha = 2/5
    ) +
    scale_x_continuous(
      limits = c(1, 7),
      breaks = c(1, 4, 7),
      minor_breaks = c(2, 3, 5, 6)
    ) +
    scale_fill_viridis_d(begin = 0, end = 0.7) +
    scale_colour_viridis_d(begin = 0, end = 0.7) +
    scale_shape_manual(
      values = c("0" = "diamond filled", "1" = "circle")
    ) +
    coord_cartesian(xlim = c(0.65, 7.35), ylim = c(0.9, 7), expand = FALSE) +
    facet_grid(. ~ ig_caste) +
    theme_grey(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(colour = "black"),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black", hjust = 0, vjust = 0),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(
      x = "Social distance"
    )
  
  # Figure 5b
  f5b <- post %>% 
    filter(outcome == "q3") %>% 
    ggplot(., aes(x = m_est, y = fct_rev(category), group = q1)) +
    stat_halfeyeh(
      aes(shape = factor(q1)),
      .width = .97, size = 0.25, fatten_point = 1,
      colour = "white", fill = "white"
    ) +
    stat_halfeyeh(
      aes(colour = ig_caste, fill = ig_caste, shape = factor(q1)),
      .width = .97, size = 0.25, fatten_point = 1,
      slab_alpha = 2/5
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 50, 100),
      minor_breaks = c(25, 75)
    ) +
    scale_fill_viridis_d(begin = 0, end = 0.7) +
    scale_colour_viridis_d(begin = 0, end = 0.7) +
    scale_shape_manual(
      values = c("0" = "diamond filled", "1" = "circle")
    ) +
    coord_cartesian(xlim = c(-5, 105), ylim = c(0.9, 7), expand = FALSE) +
    facet_grid(. ~ ig_caste) +
    theme_grey(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(colour = "black"),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black", hjust = 0, vjust = 0),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(
      x = "Feeling thermometer"
    )


# Export ------------------------------------------------------------------

  # Combine plots
  f5a + f5b + plot_layout(ncol = 1) + plot_annotation(tag_levels = "A") 
  
  # Export as .pdf
  ggsave(
    file = "figures/figure-5.pdf", 
    device = cairo_pdf,
    width = 390/.pt, height = 390/.pt, units = "mm"
  )
  
  # Export as .png
  ggsave(
    file = "figures/figure-5.png", 
    type = "cairo", dpi = 600,
    width = 390/.pt, height = 390/.pt, units = "mm"
  )
  

# Estimate ----------------------------------------------------------------

  # Summarize
  d_est <- post %>% 
    group_by_at(vars(-q1, -z_est, -m_est)) %>% 
    summarize(
      b = m_est[2] - m_est[1],
      d = z_est[2] - z_est[1]
    ) %>% 
    group_by(outcome, category, ig_caste) %>% 
    median_qi(b, d, .width = .97) %>% 
    ungroup()
  