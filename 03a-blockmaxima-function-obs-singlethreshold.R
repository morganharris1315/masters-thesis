# -------------------------------------------------------------------------
# 03a-blockmaxima-function-obs-singlethreshold.R
# -------------------------------------------------------------------------
# Feb 2026
# RX1day processing for observation data using a single threshold definition
# Threshold: value where 2/3 of annual RX1day values are above it
# -------------------------------------------------------------------------

# helper: compute RX1day with missing-data catch -------------------------
compute_RX1day_obs_st <- function(df_station, missing_day_threshold = 30) {
  # Explicit rule: if >= 30 missing daily values in a hydrological year,
  # RX1day is set to NA for that year.
  df_station %>%
    group_by(hydro_year) %>%
    summarise(
      missing_days = sum(is.na(rainfall_mm)),
      RX1day_raw = max(rainfall_mm, na.rm = TRUE),
      RX1day = ifelse(missing_days >= missing_day_threshold, NA_real_, RX1day_raw),
      .groups = "drop"
    ) %>%
    mutate(RX1day = ifelse(is.infinite(RX1day), NA_real_, RX1day)) %>%
    select(hydro_year, RX1day)
}

# Calculate single threshold ------------------------------------------------
calculate_rx1day_single_threshold_obs_st <- function(df_station_obs) {
  RX1day_df <- compute_RX1day_obs_st(df_station_obs)

  rx <- RX1day_df$RX1day
  threshold <- as.numeric(quantile(rx, probs = 1 / 3, na.rm = TRUE, type = 7))

  daily_vec <- df_station_obs$rainfall_mm
  n_days <- sum(!is.na(daily_vec))
  prop_exceed <- sum(daily_vec > threshold, na.rm = TRUE) / n_days

  list(
    threshold = threshold,
    proportion = prop_exceed
  )
}

# Single threshold legend label --------------------------------------------
make_single_threshold_label_obs_st <- function(thr) {
  c(
    "Single threshold" = paste0(
      "Single threshold (",
      round(thr$threshold, 1),
      " mm, ",
      round(thr$proportion * 100, 2),
      "%)"
    )
  )
}

# RX1day time series plot ---------------------------------------------------
plot_rx1day_obs_single_threshold_st <- function(RX1day_df, thr, labels, station_name, y_limits = NULL) {
  if (is.null(y_limits)) {
    y_max <- max(RX1day_df$RX1day, na.rm = TRUE)
    y_limits <- c(0, y_max * 1.05)
  }

  ggplot(RX1day_df, aes(x = hydro_year, y = RX1day)) +
    geom_line(color = "black") +
    geom_point(color = "black") +
    labs(
      title = paste("RX1day —", station_name),
      x = "Hydrological Year (July–June)",
      y = "RX1day (mm)",
      colour = "Threshold"
    ) +
    geom_hline(
      aes(yintercept = thr$threshold, colour = "Single threshold"),
      linetype = "solid",
      size = 1.1
    ) +
    scale_y_continuous(limits = y_limits) +
    scale_colour_manual(
      breaks = c("Single threshold"),
      labels = labels,
      values = c("Single threshold" = "#E69F00")
    ) +
    theme_thesis
}

# Histogram plot ------------------------------------------------------------
count_exceedances_per_year_obs_st <- function(df_station_obs, threshold) {
  df_station_obs %>%
    group_by(hydro_year) %>%
    summarise(
      exceedance_days = sum(rainfall_mm > threshold, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_hist_exceedances_obs_single_threshold_st <- function(df_station_obs, thr, station_name) {
  exceed_df <- count_exceedances_per_year_obs_st(df_station_obs, thr$threshold)

  x_lim <- c(-0.5, 12.5)
  y_lim <- c(0, 0.6)

  ggplot(exceed_df, aes(x = exceedance_days)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count))),
      binwidth = 1,
      boundary = 0,
      closed = "left",
      fill = "#E69F00"
    ) +
    scale_x_continuous(limits = x_lim, breaks = 0:12) +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
    theme_thesis +
    labs(
      title = station_name,
      x = "Number of exceedance days per year",
      y = "Proportion of years"
    )
}

# master function: run RX1day analysis for one station ---------------------
run_RX1day_station_analysis_obs_single_threshold_st <- function(
  combined_df_obs,
  station_name,
  output_dir = "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/RX1day_plots"
) {
  if (!station_name %in% combined_df_obs$station) {
    stop(paste("Station not found:", station_name))
  }

  df_station_obs <- combined_df_obs %>%
    filter(station == station_name) %>%
    select(observation_date, rainfall_mm, hydro_year)

  thr_station_obs <- calculate_rx1day_single_threshold_obs_st(df_station_obs)
  labels_obs <- make_single_threshold_label_obs_st(thr_station_obs)

  RX1day_df <- compute_RX1day_obs_st(df_station_obs)

  p_ts <- plot_rx1day_obs_single_threshold_st(
    RX1day_df = RX1day_df,
    thr = thr_station_obs,
    labels = labels_obs,
    station_name = station_name
  )

  p_hist <- plot_hist_exceedances_obs_single_threshold_st(
    df_station_obs = df_station_obs,
    thr = thr_station_obs,
    station_name = station_name
  )

  station_safe <- gsub(" ", "_", station_name)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Time series
  ggsave(
    filename = file.path(output_dir, paste0(station_safe, "_RX1day_TimeSeries.png")),
    plot = p_ts,
    width = fig_width_full,
    height = fig_height_tall,
    dpi = 300
  )

  # Histogram
  ggsave(
    filename = file.path(output_dir, paste0(station_safe, "_RX1day_Histogram.png")),
    plot = p_hist,
    width = fig_width_hoz_full,
    height = fig_height_med,
    dpi = 300
  )

  invisible(list(
    timeseries = p_ts,
    histogram = p_hist
  ))
}

# Usage --------------------------------------------------------------------

# all regions; save into region/single_threshold subfolders
region_output_dirs_obs_single_threshold_st <- c(
  coromandel = glue("{base_raw_dir}/obs_data/coromandel/single_threshold"),
  far_north = glue("{base_raw_dir}/obs_data/far_north/single_threshold"),
  top_of_south = glue("{base_raw_dir}/obs_data/top_of_south/single_threshold"),
  waikato = glue("{base_raw_dir}/obs_data/waikato/single_threshold")
)

imap(region_output_dirs_obs_single_threshold_st, function(region_dir, region_name) {
  stations <- unique(combined_df_obs$station[combined_df_obs$region == region_name])

  walk(
    stations,
    ~ run_RX1day_station_analysis_obs_single_threshold_st(
      combined_df_obs,
      station_name = .x,
      output_dir = region_dir
    )
  )
})

head(combined_df_obs)
