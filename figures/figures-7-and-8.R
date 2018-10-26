rm(list = ls())

# Library -----------------------------------------------------------------
  library(tidyverse); library(tidybayes)


# Import ------------------------------------------------------------------
  post <- read_rds("models/ldps_post.rds") 
  
  # Prepare for figures
  post <- post %>%
    separate(item, c("item", "group"), sep = "_") %>%
    mutate(
      group = recode_factor(
        group,
        "ig" = "...from your\nown background",
        "scst" = "SC/ST",
        "obc" = "OBC",
        "gm" = "GM",
        "hindu" = "Hindus",
        "muslim" = "Muslims",
        .ordered = TRUE
      ),
      ig_caste = recode_factor(
        ig_caste,
        "scst" = "SC/ST",
        "obc" = "OBC",
        "gm" = "GM",
        .ordered = TRUE
      )
    )
  

# Figure 7 ----------------------------------------------------------------
  # Labels for Figure 7
  d_lbl <- tibble(
    x = c(2.5, 5.5), xend = c(1.5, 6.5),
    ig_caste = "SC/ST", group = ordered("Muslims", levels = unique(post$group)),
    label = c("harder", "easier")
  )
  
  # Figure 7
  post %>%
    filter(item == "ld") %>%
  ggplot(., aes(x = m_est, y = ig_caste)) +
    geom_vline(xintercept = 4, colour = "white", size = 1) +
    geom_halfeyeh(
      aes(shape = "estimate"),
      colour = "white", fill = "white",
      point_interval = median_hdi, .width = 0.97,
      rel_min_height = 0.01, size = 0
    ) +
    geom_halfeyeh(
      aes(colour = ig_caste, fill = ig_caste, shape = "estimate"),
      point_interval = median_hdi, .width = 0.97,
      rel_min_height = 0.01, size = 0.5, fatten_point = 2.5
    ) +
    geom_segment(
      data = d_lbl, aes(x = x, xend = xend, yend = ig_caste),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed")
    ) +
    geom_text(
      data = d_lbl, aes(label = label), x = c(2, 6), y = c(1, 1),
      size = 10*0.8/.pt, vjust = -0.5
    ) +
    scale_x_continuous(breaks = 1:7, minor_breaks = NULL) +
    scale_fill_viridis_d(begin = 0, end = 0.7, direction = -1, alpha = 2/5) +
    scale_colour_viridis_d(begin = 0, end = 0.7, direction = -1) +
    scale_shape_manual(values = c(18)) +
    coord_fixed(ratio = 1/2, xlim = c(0.9,7.1), ylim = c(0.8, 4), expand = FALSE) +
    facet_grid(group ~ .) +
    labs(
      x = "Perceived (dis-)advantage",
      y = "Ingroup"
    ) +
    theme_grey(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(vjust = 0),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black", hjust = 0, vjust = 0),
      axis.ticks.y = element_blank()
    )
  
  # Export
  ggsave(file = "figures/figure-7.pdf", width = 390/.pt, height = 160, units = "mm")

# Figure 8 ----------------------------------------------------------------
  # Labels for Figure 8
  d_lbl <- tibble(
    x = c(2.25, 3.75), xend = c(1.25, 4.75),
    ig_caste = "SC/ST", group = ordered("Muslims", levels = unique(post$group)),
    label = c("oppose", "support")
  )
  
  # Figure 8
  post %>%
    filter(item == "ps") %>%
    ggplot(., aes(x = m_est, y = ig_caste)) +
    geom_vline(xintercept = 3, colour = "white", size = 1) +
    geom_halfeyeh(
      aes(shape = "estimate"),
      colour = "white", fill = "white",
      point_interval = median_hdi, .width = 0.97,
      rel_min_height = 0.01, size = 0
    ) +
    geom_halfeyeh(
      aes(colour = ig_caste, fill = ig_caste, shape = "estimate"),
      point_interval = median_hdi, .width = 0.97,
      rel_min_height = 0.01, size = 0.5, fatten_point = 2.5
    ) +
    geom_segment(
      data = d_lbl, aes(x = x, xend = xend, yend = ig_caste),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed")
    ) +
    geom_text(
      data = d_lbl, aes(label = label), x = c(1.75, 4.25), y = c(1, 1),
      size = 10*0.8/.pt, vjust = -0.5
    ) +
    scale_x_continuous(breaks = 1:5, minor_breaks = NULL) +
    scale_fill_viridis_d(begin = 0, end = 0.7, direction = -1, alpha = 2/5) +
    scale_colour_viridis_d(begin = 0, end = 0.7, direction = -1) +
    scale_shape_manual(values = c(18)) +
    coord_fixed(ratio = 1/3, xlim = c(0.9, 5.1), ylim = c(0.8, 4), expand = FALSE) +
    facet_grid(group ~ .) +
    labs(
      x = "Policy support",
      y = "Ingroup"
    ) +
    theme_grey(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(vjust = 0),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_text(colour = "black", hjust = 0, vjust = 0),
      axis.ticks.y = element_blank()
    )
  
  # Export
  ggsave(file = "figures/figure-8.pdf", width = 390/.pt, height = 80, units = "mm")
   