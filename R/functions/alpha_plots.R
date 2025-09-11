alpha_plots <- function(
  df,
  .x,
  .y,
  .color = .x,
  .facet = NULL,
  title,
  subtitle,
  x_lab,
  y_lab,
  legend_by,
  alpha
) {
  # Capture quosures safely
  .x <- rlang::enquo(.x)
  .y <- rlang::enquo(.y)
  .facet <- rlang::enquo(.facet)
  .color <- rlang::enquo(.color)

  plot <- ggplot(
    df,
    aes(x = !!.x, y = !!.y)
  ) +
    geom_jitter(aes(color = !!.color), width = 0.2, alpha = alpha) +
    geom_boxplot(alpha = 0, width = 0.6) +

    if (!rlang::quo_is_null(.facet)) {
      plot <- plot +
        facet_wrap(vars(!!.facet), scales = "free_y")
    }
  plot <- plot +
    theme_bw() +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_lab,
      y = y_lab,
      color = legend_by
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

  return(plot)
}
