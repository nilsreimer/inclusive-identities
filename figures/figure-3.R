rm(list = ls())

# Library -----------------------------------------------------------------

  # Load packages
  library(tidyverse); library(patchwork); library(tidybayes)


# Import ------------------------------------------------------------------
  
  # Observed proportions
  dl <- read_rds("data/d1.rds") %>% 
    filter(!is.na(q1)) %>%
    mutate(
      ig_caste = recode_factor(
        ig_caste, 
        "gm" = "GM", 
        "obc" = "OBC",
        "scst" = "SC/ST"
      ),
      caste = recode_factor(
        caste, 
        "GM" = "GM", 
        "OBC" = "OBC", 
        "SCST" = "SC/ST", 
        "N/A" = NA_character_
      )
    )
  
  # Estimated probabilites
  post <- read_rds("results/q1_post.rds") %>%
    unnest(p_est) %>%
    mutate(
      ig_caste = recode_factor(
        ig_caste,
        "1" = "GM",
        "2" = "OBC",
        "3" = "SC/ST"
      )
    )
  

# Summarise ---------------------------------------------------------------
  
  # Observed proportions
  d_obs <- dl %>% 
    group_by(target, category, lname, nation, religion, caste, ig_caste) %>%
    summarise(p_obs = mean(q1))
  
  # Estimated probabilities
  d_est <- post %>%
    rename(p_est = p_k) %>%
    group_by(model, category, ig_caste) %>%
    median_qi(.width = c(.67, .89, .97))


# Plot --------------------------------------------------------------------
  
  # Make plot
  p3 <- d_est %>%
    left_join(d_obs, ., by = c("category", "ig_caste")) %>%
    filter(model == 3) %>%
    mutate(column = ig_caste) %>% 
    # mutate(
    #   column = factor(
    #     ig_caste,
    #     levels = c("GM", "OBC", "SC/ST", "Combined")
    #   )
    # ) %>% 
    # bind_rows(
    #   d_est %>%
    #     left_join(d_obs, ., by = c("category", "ig_caste")) %>%
    #     filter(model == 3) %>% 
    #     mutate(
    #       column = factor(
    #         "Combined",
    #         levels = c("GM", "OBC", "SC/ST", "Combined")
    #       ),
    #       p_obs = NA_real_
    #     )
    # ) %>% 
  ggplot(., aes(x = p_est, y = 25 - target, xmin = .lower, xmax = .upper)) +
    annotate(
      geom = "rect",
      xmin = c(-Inf, 1),
      xmax = c(0, Inf),
      ymin = c(-Inf, -Inf),
      ymax = c(Inf, Inf),
      colour = NA,
      fill = "white"
    ) +
    geom_ribbon(
      aes(group = interaction(1-.width, ig_caste)),
      fill = "white"
    ) +
    geom_ribbon(
      aes(fill = ig_caste, group = interaction(1-.width, ig_caste)),
      alpha = 1/3
    ) +
    geom_path(aes(colour = ig_caste), size = 0.25) +
    geom_point(aes(x = p_obs), shape = 3, colour = "grey60", size = 0.75) +
    geom_point(aes(colour = ig_caste), size = 0.75) +
    geom_hline(
      yintercept = c(4.5, 8.5, 12.5, 16.5, 20.5), 
      colour = "white", 
      size = 2
    ) +
    geom_hline(yintercept = c(0.5, 24.5), colour = "black") +
    scale_x_continuous(
      breaks = seq(0.0, 1.0, 0.5), 
      minor_breaks = NULL
    ) +
    scale_y_continuous(
      breaks = seq(1, 24, 1), 
      minor_breaks = NULL
    ) +
    scale_colour_viridis_d(begin = 0.0, end = 0.7) +
    scale_fill_viridis_d(begin = 0.0, end = 0.7) +
    coord_cartesian(
      xlim = c(-0.1, 1.1), 
      ylim = c(0.5, 24.5), 
      expand = FALSE
    ) +
    facet_grid(. ~ column) +
    theme_grey(base_size = 10, base_line_size = 0.25) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(colour = "black", size = 10),
      axis.text.x = element_text(colour = "black", size = 9),
      axis.text.y = element_blank(),
      axis.title.x = element_text(colour = "black", size = 9),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(colour = "black"),
      axis.ticks.y = element_blank()
    ) +
    labs(
      x = 'Pr("us"|M3)'
    )


# Table -------------------------------------------------------------------
  
  # Make table
  t3 <- dl %>%
    distinct(target, nation, religion, caste) %>%
    arrange(target) %>% 
    mutate(
      `#` = scales::number(target, accuracy = 1, trim = FALSE),
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
    pivot_longer(
      -target,
      names_to = "column",
      values_to = "text"
    ) %>% 
    mutate(
      column = recode_factor(
        column,
        "nation" = "Nation",
        "religion" = "Religion",
        "caste" = "Caste",
        "#" = "#"
      ),
      .lower = if_else(column == "#", -0.025, -0.01),
      .upper = if_else(column == "#",  0.025,  0.25)
    ) %>% 
  ggplot(., aes(x = 0, y = 25 - target, label = text)) +
    geom_linerange(aes(xmin = .lower, xmax = .upper), colour = NA) +
    geom_text(aes(hjust = rep(c(0, 0, 0, 0.5), 24)), size = 10/.pt) +
    geom_hline(
      yintercept = c(0.5, 24.5),
      colour = "black"
    ) +
    geom_hline(
      yintercept = c(4.5, 8.5, 12.5, 16.5, 20.5), 
      colour = "black", 
      size = 0.125
    ) +
    facet_grid(. ~ column, scales = "free_x", space = "free_x") +
    scale_y_continuous(limits = c(0.5, 24.5), expand = c(0, 0)) +
    theme_grey(base_size = 10, base_line_size = 0.25) +
    theme(
      strip.background = element_rect(fill = NA),
      strip.text = element_text(colour = "black", size = 10, hjust = 0),
      panel.spacing = unit(0.0, "lines"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )


# Export ------------------------------------------------------------------

  # Combine table and plot  
  t3 + p3 + plot_layout(widths = c(0.52, 0.48))
    
  # Export as .pdf
  ggsave(
    file = "figures/figure-3.pdf", 
    device = cairo_pdf,
    width = 390/.pt, height = 160, units = "mm"
  )
  
  # Export as .png
  ggsave(
    file = "figures/figure-3.png", 
    type = "cairo", dpi = 600,
    width = 390/.pt, height = 160, units = "mm"
  )


# Estimate ----------------------------------------------------------------

  # Extract estimates for main text
  post %>% 
    filter(model == 3L) %>% 
    group_by(category, ig_caste) %>% 
    median_qi(p_k, .width = 0.97) %>% 
    ungroup() %>% 
    mutate_at(
      vars(p_k, .lower, .upper),
      ~round(., digits = 2) %>% 
        scales::number(accuracy = 0.01) %>% 
        str_replace("0", "")
    ) %>% 
    mutate(
      text = glue::glue("Pr("us"|M3) = {p_k}, [{.lower}, {.upper}]")
    )
  
  # Extract differences for main text
  post %>% 
    filter(model == 3L) %>% 
    group_by(ig_caste) %>% 
    pivot_wider(
      names_from = category,
      names_prefix = "p_",
      values_from = p_k
    ) %>% 
    mutate(
      d = p_2 - p_3,
      r = p_2 / p_3
    ) %>% 
    median_qi(d, r, .width = 0.97) %>% 
    mutate_if(is.numeric, round, digits = 2)
  
  