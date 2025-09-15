acf_preprocess <- function(
  df,
  feature,
  time_col = "sampling_date",
  agg_fun = mean,
  na_rm = TRUE
) {
  series <- df %>%
    filter(!is.na(.data[[time_col]])) %>%
    group_by(.data[[time_col]]) %>%
    summarise(
      value = agg_fun(
        as.numeric(.data[[feature]]),
        na.rm = na_rm
      ),
      .groups = "drop"
    ) %>%
    arrange(.data[[time_col]]) %>%
    filter(is.finite(value))

  return(series)
}


acf_plots <- function(
  df,
  feature,
  time_col = "sampling_date",
  title = NULL,
  lag_max = NULL,
  agg_fun = mean
) {
  series <- acf_preprocess(
    df,
    feature = feature,
    time_col = time_col,
    agg_fun = agg_fun
  )
  ac <- stats::acf(series$value, plot = FALSE, na.action = stats::na.pass)

  df_ac <- data.frame(
    lag = as.integer(ac$lag[, 1, 1]),
    acf = as.numeric(ac$acf[, 1, 1])
  )

  conf <- 1.96 / sqrt(ac$n.used)
  if (is.null(title)) {
    title <- "ACF"
  }

  ggplot(df_ac, aes(lag, acf)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_hline(
      yintercept = c(-conf, conf),
      linetype = "dashed",
      color = "red"
    ) +
    geom_col(width = 0.8, fill = "steelblue") +
    labs(title = title, x = "Lag", y = "ACF") +
    theme_minimal()
}
#----------------
