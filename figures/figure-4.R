rm(list = ls())

# Library -----------------------------------------------------------------
    
  # Load packages
  library(tidyverse); library(tidybayes)


# Import ------------------------------------------------------------------
  
  # Estimated probabilites
  post <- read_rds("results/q1_ic_post.rds")
  

# Summarise ---------------------------------------------------------------
  
  # Estimated probabilities
  d_est <- post %>%
    group_by(category, ig_caste, item, response) %>%
    median_hdi(p_est, .width = c(.67, .89, .97)) %>%
    ungroup() %>%
    mutate(
      ig_caste = case_when(
        ig_caste == 1 ~ "GM",
        ig_caste == 2 ~ "OBC",
        ig_caste == 3 ~ "SC/ST"
      ),
      category = recode_factor(
        category,
        "1" = "Hindu, GM",
        "2" = "Hindu, OBC",
        "3" = "Hindu, SC/ST",
        "4" = "Muslim, OBC",
        .ordered = TRUE
      ),
      item = case_when(
        item == "nc" ~ "Negative contact",
        item == "of" ~ "Outgroup friendship"
      )
    ) %>% 
    distinct()


# Figure ------------------------------------------------------------------

  # Make figure
  d_est %>%
    ggplot(., aes(
      x = response, y = p_est, ymin = .lower, ymax = .upper, 
      colour = ig_caste, fill = ig_caste
    )) +
    geom_ribbon(aes(group = interaction(.width, ig_caste)), colour = NA, fill = "white") +
    geom_ribbon(aes(group = interaction(.width, ig_caste)), colour = NA, alpha = 0.20) +
    geom_line() +
    geom_point(aes(shape = ig_caste), size = 0.75) +
    scale_x_continuous(limits = c(1, 5), breaks = 1:5, minor_breaks = 0) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
    scale_colour_viridis_d(begin = 0.0, end = 0.7) +
    scale_fill_viridis_d(begin = 0.0, end = 0.7) +
    scale_shape_manual(values = c("circle", "triangle", "square", "diamond filled")) +
    labs(
      x = "Amount",
      y ='Pr("us"|M5)',
      colour = "Ingroup:",
      fill = "Ingroup:",
      shape = "Ingroup:"
    ) +
    facet_grid(item ~ category) +
    coord_fixed(4) +
    theme_grey(base_size = 10, base_line_size = 0.25) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = NA, colour = NA),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black")
    )
  

# Export ------------------------------------------------------------------

  # Export as .pdf
  ggsave(
    file = "figures/figure-4.pdf", 
    device = cairo_pdf,
    width = 390/.pt, height = 90, units = "mm"
  )
  
  # Export as .png
  ggsave(
    file = "figures/figure-4.png", 
    type = "cairo", dpi = 600,
    width = 390/.pt, height = 90, units = "mm"
  )
  