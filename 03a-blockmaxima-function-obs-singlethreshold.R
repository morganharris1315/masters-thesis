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
      values = c("Single threshold" = "#93acff")
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
      fill = "#93acff"
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


# Example years ------------------------------------------------------------

# Select example years safely
select_example_years_obs_st <- function(df_station, threshold) {
  # Compute RX1day per year
  RX1day_df <- compute_RX1day_obs_st(df_station) %>%
    rename(Year = hydro_year)

  # Count exceedances per year
  exceedances <- df_station %>%
    group_by(hydro_year) %>%
    summarise(exceedance_days = sum(rainfall_mm > threshold, na.rm = TRUE), .groups = "drop") %>%
    rename(Year = hydro_year)

  summary_df <- RX1day_df %>%
    left_join(exceedances, by = "Year") %>%
    filter(!is.na(RX1day))

  muted_year <- summary_df %>% filter(exceedance_days == 0) %>% arrange(RX1day) %>% slice(1) %>% pull(Year)
  if (length(muted_year) == 0) muted_year <- NA_integer_

  single_year <- summary_df %>% filter(exceedance_days == 1) %>% arrange(desc(RX1day)) %>% slice(1) %>% pull(Year)
  if (length(single_year) == 0) single_year <- NA_integer_

  high_year <- summary_df %>% arrange(desc(exceedance_days), desc(RX1day)) %>% slice(1) %>% pull(Year)
  if (length(high_year) == 0) high_year <- NA_integer_

  list(
    muted = muted_year,
    single = single_year,
    high = high_year
  )
}

# Plot daily timeseries for one year
plot_daily_example_obs_st <- function(df_station, year, threshold, title, y_max, colour, y_lab = "Daily Rainfall (mm)") {
  if (is.na(year)) return(NULL)

  ts_df <- df_station %>%
    filter(hydro_year == year)

  ggplot(ts_df, aes(x = observation_date, y = rainfall_mm, group = 1)) +
    geom_line() +
    geom_point(size = 1) +
    geom_hline(
      yintercept = threshold,
      linetype = "dashed",
      colour = colour,
      size = 1.1
    ) +
    scale_x_date(
      date_labels = "%b",
      date_breaks = "3 months",
      expand = expansion(add = c(0, 0))
    ) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = 50)) +
    labs(title = title, x = "Date", y = y_lab) +
    theme_thesis +
    theme(plot.title = element_text(hjust = 0.5))
}

# Plot example years for one station
plot_exceedance_examples_obs_st <- function(df_station, station_name, threshold, colour = "#93acff") {
  years <- select_example_years_obs_st(df_station, threshold)

  # Determine max y across valid years
  valid_years <- years[!is.na(unlist(years))]
  if (length(valid_years) == 0) {
    warning(paste("No valid example years for station:", station_name))
    return(NULL)
  }

  y_max <- df_station %>%
    filter(hydro_year %in% unlist(valid_years)) %>%
    summarise(max_val = max(rainfall_mm, na.rm = TRUE)) %>%
    pull(max_val) %>% {
      . + 10
    }

  # Create plots safely
  plots <- list()

  if (!is.na(years$muted)) {
    plots$muted <- plot_daily_example_obs_st(df_station, years$muted, threshold, paste0("Muted Year (", years$muted, ")"), y_max, colour)
  }
  if (!is.na(years$single)) {
    plots$single <- plot_daily_example_obs_st(df_station, years$single, threshold, paste0("Single Exceedance Year (", years$single, ")"), y_max, colour)
  }
  if (!is.na(years$high)) {
    plots$high <- plot_daily_example_obs_st(df_station, years$high, threshold, paste0("High Exceedance Year (", years$high, ")"), y_max, colour)
  }

  # Combine available plots
  wrap_plots(plots) +
    plot_annotation(title = station_name, theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
}

# Master function for example years per station
run_example_years_obs_single_threshold_st <- function(df_obs, station_name, output_dir) {
  if (!station_name %in% df_obs$station) {
    stop(paste("Station not found:", station_name))
  }

  df_station <- df_obs %>% filter(station == station_name)

  # Calculate single threshold
  thr_station <- calculate_rx1day_single_threshold_obs_st(df_station)

  # Example-year plot (single threshold)
  p_example <- plot_exceedance_examples_obs_st(df_station, station_name, thr_station$threshold)

  # Save plot if it exists
  station_safe <- gsub(" ", "_", station_name)

  if (!is.null(p_example)) {
    ggsave(
      filename = file.path(output_dir, paste0(station_safe, "_RX1day_ExampleYears_SingleThreshold.png")),
      plot = p_example,
      width = fig_width_full,
      height = fig_height_tall,
      dpi = 300
    )
  }

  invisible(list(single_threshold = p_example))
}

imap(region_output_dirs_obs_single_threshold_st, function(region_dir, region_name) {
  # Get all stations in the region
  stations <- unique(combined_df_obs$station[combined_df_obs$region == region_name])

  # Loop through stations
  walk(stations, function(station_name) {
    # Subset data for this station
    df_station <- combined_df_obs %>% filter(station == station_name)

    # Calculate threshold and check if valid example years exist
    thr_station <- calculate_rx1day_single_threshold_obs_st(df_station)
    years_single <- select_example_years_obs_st(df_station, thr_station$threshold)
    all_years <- unlist(years_single)

    # Skip station if no valid years
    if (all(is.na(all_years))) {
      message(paste("Skipping station (no valid example years):", station_name))
      return(NULL)
    }

    # Otherwise, run example-years analysis
    run_example_years_obs_single_threshold_st(combined_df_obs, station_name, region_dir)
  })
})


# Exceedance summary table -------------------------------------------------
build_high_exceedance_table_obs_single_threshold_st <- function(combined_df_obs) {
  stations <- unique(combined_df_obs$station)

  bind_rows(
    lapply(stations, function(station_name) {
      df_station <- combined_df_obs %>% filter(station == station_name)
      region <- unique(df_station$region)

      # Calculate station-specific single threshold
      thr <- calculate_rx1day_single_threshold_obs_st(df_station)$threshold

      # RX1day per year
      RX1day_df <- compute_RX1day_obs_st(df_station) %>%
        rename(Year = hydro_year)

      # Exceedances per year
      exc_df <- count_exceedances_per_year_obs_st(df_station, thr) %>%
        rename(
          Exceedance_Days = exceedance_days,
          Year = hydro_year
        )

      # Annual rainfall per year
      annual_rain_df <- df_station %>%
        group_by(hydro_year) %>%
        summarise(AnnualRain = sum(rainfall_mm, na.rm = TRUE), .groups = "drop") %>%
        rename(Year = hydro_year)

      RX1day_df %>%
        left_join(exc_df, by = "Year") %>%
        left_join(annual_rain_df, by = "Year") %>%
        filter(!is.na(Exceedance_Days)) %>%
        arrange(desc(RX1day)) %>%
        mutate(
          RX1day_Rank = row_number(),
          RX1day_Percentile = (1 - (RX1day_Rank - 1) / (n() - 1)) * 100,
          Threshold = "Single threshold",
          Station = station_name,
          Region = region
        ) %>%
        select(
          Region,
          Station,
          Threshold,
          Year,
          Exceedance_Days,
          RX1day,
          RX1day_Rank,
          RX1day_Percentile,
          AnnualRain
        )
    })
  ) %>%
    arrange(desc(Exceedance_Days), desc(RX1day))
}

high_exceedance_table_single_threshold <-
  build_high_exceedance_table_obs_single_threshold_st(combined_df_obs)

# Top-20 exceedance years across all stations (no filtering required)
top20_single_threshold <- high_exceedance_table_single_threshold %>%
  slice_max(order_by = Exceedance_Days, n = 20, with_ties = FALSE)

print(top20_single_threshold, n = 20)

obs_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/obs_data"

write.csv(
  top20_single_threshold,
  file = file.path(obs_dir, "highexceedance_singlethreshold_obs.csv"),
  row.names = FALSE
)
