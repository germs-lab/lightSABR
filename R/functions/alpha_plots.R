alpha_plots <- function(
  df,
  feat,
  .x,
  .y,
  .color = .x,
  .facet,
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
  feat <- rlang::enquo(feat)
  .facet <- rlang::enquo(.facet)
  .color <- rlang::enquo(.color)

  plot <- ggplot(
    df,
    aes(x = !!.x, y = !!.y)
  ) +
    geom_jitter(aes(color = !!.color), width = 0.2, alpha = alpha) +
    geom_boxplot(alpha = 0, width = 0.6) +
    facet_wrap(vars(!!.facet), scales = "free_y") +
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
