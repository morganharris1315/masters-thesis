# -------------------------------------------------------------------------
# 03a-blockmaxima-function-obs-singlethreshold.R
# -------------------------------------------------------------------------
# Jan 2026
# Creating functions to make the processing of observation data easier
# -------------------------------------------------------------------------

# function: compute RX1day metrics & thresholds -------------------------
calculate_rx1day_thresholds_obs <- function(df_station_obs) {
  
  RX1day_df <- df_station_obs %>%
    group_by(hydro_year) %>%
    summarise(
      RX1day = max(rainfall_mm, na.rm = TRUE),
      .groups = "drop"
    )
  
  rx <- RX1day_df$RX1day
  
  min_rx  <- min(rx, na.rm = TRUE)
  mean_rx <- mean(rx, na.rm = TRUE)
  max_rx  <- max(rx, na.rm = TRUE)
  mid_rx  <- (min_rx + mean_rx) / 2
  
  daily_vec <- df_station_obs$rainfall_mm
  n_days <- sum(!is.na(daily_vec))
  
  prop_exceed <- c(
    min  = sum(daily_vec > min_rx,  na.rm = TRUE) / n_days,
    mid  = sum(daily_vec > mid_rx,  na.rm = TRUE) / n_days,
    mean = sum(daily_vec > mean_rx, na.rm = TRUE) / n_days,
    max  = sum(daily_vec > max_rx,  na.rm = TRUE) / n_days
  )
  
  prop_onein360 <- 1 / 360
  
  threshold_onein360 <- quantile(
    daily_vec,
    probs = 1 - prop_onein360,
    na.rm = TRUE,
    type = 7
  )
  
  list(
    thresholds = c(
      min      = min_rx,
      mid      = mid_rx,
      mean     = mean_rx,
      max      = max_rx,
      onein360 = as.numeric(threshold_onein360)
    ),
    proportions = c(
      prop_exceed,
      onein360 = prop_onein360
    )
  )
}

#  function: RX1day time series plot --------------------------------------
make_labels_obs <- function(thr) {
  labels <- c("Min RX1day" = paste0("Min RX1day (", round(thr$thresholds["min"], 1), " mm, ", round(thr$proportions["min"] * 100, 2), "%)"),
              "Mid RX1day" = paste0("Mid RX1day (", round(thr$thresholds["mid"], 1), " mm, ", round(thr$proportions["mid"] * 100, 2), "%)"),
              "Mean RX1day" = paste0("Mean RX1day (", round(thr$thresholds["mean"], 1), " mm, ", round(thr$proportions["mean"] * 100, 2), "%)"),
              "1 in 360 day threshold" = paste0("1 in 360 day (", round(thr$thresholds["onein360"], 1), " mm, ", round(thr$proportions["onein360"] * 100, 2), "%)"),
              "Max RX1day" = paste0("Max RX1day (", round(thr$thresholds["max"], 1), " mm, ",round(thr$proportions["max"] * 100, 2), "%)"))
  return(labels)
}

plot_rx1day_obs <- function(RX1day_df, thr, labels, station_name, y_limits = NULL) {
  if (is.null(y_limits)) {y_max <- max(RX1day_df$RX1day, na.rm = TRUE) 
  y_limits <- c(0, y_max * 1.05)}
  ggplot(RX1day_df, aes(x = hydro_year, y = RX1day)) + geom_line(color = "black") + geom_point(color = "black") +
    labs(title = paste("RX1day —", station_name), x = "Hydrological Year (July–June)", y = "RX1day (mm)", colour = "Threshold") +
    geom_hline(aes(yintercept = thr$thresholds["min"], colour = "Min RX1day"), linetype = "dashed", size = 1.1) +
    geom_hline(aes(yintercept = thr$thresholds["mid"], colour = "Mid RX1day"), linetype = "dashed", size = 1.1) +
    geom_hline(aes(yintercept = thr$thresholds["mean"], colour = "Mean RX1day"), linetype = "dashed", size = 1.1) +
    geom_hline(aes(yintercept = thr$thresholds["max"], colour = "Max RX1day"), linetype = "dashed", size = 1.1) +
    geom_hline(aes(yintercept = thr$thresholds["onein360"], colour = "1 in 360 day threshold"), linetype = "solid", size = 1.1) +
    scale_y_continuous(limits = y_limits) +
    scale_colour_manual(
      breaks = c("Max RX1day", "Mean RX1day", "1 in 360 day threshold", "Mid RX1day", "Min RX1day"),
      labels = labels,
      values = c(
        "Min RX1day"             = "#56B4E9",
        "Mid RX1day"             = "#0072B2",
        "1 in 360 day threshold" = "#E69F00",
        "Mean RX1day"            = "#009E73",
        "Max RX1day"             = "#D55E00"
      )
    ) +
    theme_thesis
}

#  functions: histograms --------------------------------------
count_exceedances_per_year_obs <- function(df_station_obs, threshold) {
  
  df_station_obs %>%
    group_by(hydro_year) %>%
    summarise(
      exceedance_days = sum(rainfall_mm > threshold, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_hist_exceedances_obs <- function(df_station_obs, thr, station_name) {
  
  mid_df  <- count_exceedances_per_year_obs(df_station_obs, thr$thresholds["mid"])
  mean_df <- count_exceedances_per_year_obs(df_station_obs, thr$thresholds["mean"])
  d360_df <- count_exceedances_per_year_obs(df_station_obs, thr$thresholds["onein360"])
  
  x_lim <- c(-0.5, 12.5)
  y_lim <- c(0, 0.6)
  
  p_mid <- ggplot(mid_df, aes(x = exceedance_days)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count))),
      binwidth = 1, boundary = 0, closed = "left",
      fill = "#0072B2"
    ) +
    scale_x_continuous(limits = x_lim, breaks = 0:12) +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
    theme_thesis +
    labs(
      title = "> Mid RX1day",
      x = "Number of exceedance days per year",
      y = "Proportion of years"
    )
  
  p_360 <- ggplot(d360_df, aes(x = exceedance_days)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count))),
      binwidth = 1, boundary = 0, closed = "left",
      fill = "#E69F00"
    ) +
    scale_x_continuous(limits = x_lim, breaks = 0:12) +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
    theme_thesis +
    labs(
      title = "> 1-in-360 day threshold",
      x = "Number of exceedance days per year",
      y = "Proportion of years"
    )
  
  p_mean <- ggplot(mean_df, aes(x = exceedance_days)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count))),
      binwidth = 1, boundary = 0, closed = "left",
      fill = "#009E73"
    ) +
    scale_x_continuous(limits = x_lim, breaks = 0:12) +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
    theme_thesis + 
    labs(
      title = "> Mean RX1day",
      x = "Number of exceedance days per year",
      y = "Proportion of years"
    )
  
  (p_mid | p_360 | p_mean) +
    plot_annotation(
      title = station_name,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    )
}


# master function: run RX1day analysis for one station --------------------
run_RX1day_station_analysis_obs <- function(combined_df_obs,station_name,
                                            output_dir = "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/RX1day_plots") {
  
  if (!station_name %in% combined_df_obs$station) {
    stop(paste("Station not found:", station_name))
  }
  
  df_station_obs <- combined_df_obs %>%
    filter(station == station_name) %>%
    select(observation_date, rainfall_mm, hydro_year)
  
  thr_station_obs <- calculate_rx1day_thresholds_obs(df_station_obs)
  labels_obs <- make_labels_obs(thr_station_obs)
  
  RX1day_df <- df_station_obs %>%
    group_by(hydro_year) %>%
    summarise(RX1day = max(rainfall_mm, na.rm = TRUE), .groups = "drop")
  
  p_ts <- plot_rx1day_obs(
    RX1day_df = RX1day_df,
    thr = thr_station_obs,
    labels = labels_obs,
    station_name = station_name
  )
  
  p_hist <- plot_hist_exceedances_obs(
    df_station_obs = df_station_obs,
    thr = thr_station_obs,
    station_name = station_name
  )
  
  station_safe <- gsub(" ", "_", station_name)
  
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
    histogram  = p_hist
  ))
}


# Usage -------------------------------------------------------------------

#all reigions
region_output_dirs <- c(
  coromandel   = glue("{base_raw_dir}/obs_data/coromandel"),
  far_north    = glue("{base_raw_dir}/obs_data/far_north"),
  top_of_south = glue("{base_raw_dir}/obs_data/top_of_south"),
  waikato      = glue("{base_raw_dir}/obs_data/waikato")
)

imap(region_output_dirs, function(region_dir, region_name) {
  
  # Get all stations in this region
  stations <- unique(combined_df_obs$station[combined_df_obs$region == region_name])
  
  # Run RX1day analysis for each station in its region folder
  walk(
    stations,
    ~ run_RX1day_station_analysis_obs(
      combined_df_obs,
      station_name = .x,
      output_dir = region_dir 
    )
  )
})


head(combined_df_obs)


# Example years -----------------------------------------------------------

# Compute RX1day per year
compute_RX1day_obs <- function(df_station) {
  df_station %>%
    group_by(hydro_year) %>%
    summarise(RX1day = max(rainfall_mm, na.rm = TRUE), .groups = "drop") %>%
    rename(Year = hydro_year)
}

# Select example years safely
select_example_years_obs <- function(df_station, threshold) {
  
  # Compute RX1day per year
  RX1day_df <- compute_RX1day_obs(df_station)
  
  # Count exceedances per year
  exceedances <- df_station %>%
    group_by(hydro_year) %>%
    summarise(exceedance_days = sum(rainfall_mm > threshold, na.rm = TRUE), .groups = "drop") %>%
    rename(Year = hydro_year)
  
  summary_df <- RX1day_df %>%
    left_join(exceedances, by = "Year")
  
  muted_year  <- summary_df %>% filter(exceedance_days == 0) %>% arrange(RX1day) %>% slice(1) %>% pull(Year)
  if(length(muted_year) == 0) muted_year <- NA_integer_
  
  single_year <- summary_df %>% filter(exceedance_days == 1) %>% arrange(desc(RX1day)) %>% slice(1) %>% pull(Year)
  if(length(single_year) == 0) single_year <- NA_integer_
  
  high_year   <- summary_df %>% arrange(desc(exceedance_days), desc(RX1day)) %>% slice(1) %>% pull(Year)
  if(length(high_year) == 0) high_year <- NA_integer_
  
  list(
    muted  = muted_year,
    single = single_year,
    high   = high_year
  )
}

# Extract daily timeseries for a specific year
extract_daily_timeseries_obs <- function(df_station, year) {
  df_station %>%
    filter(hydro_year == year) %>%
    mutate(day = as.integer(format(observation_date, "%j"))) %>%  # day of year
    select(day, rainfall_mm)
}

# Plot daily timeseries for one year
plot_daily_example_obs <- function(df_station, year, threshold, title, y_max, colour, y_lab = "Daily Rainfall (mm)") {
  
  if(is.na(year)) return(NULL)  # skip missing years safely
  
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
    )+
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = 50)) +
    labs(title = title, x = "Date", y = y_lab) +
    theme_thesis + 
    theme(plot.title = element_text(hjust = 0.5))
}


# Plot example years for one station
plot_exceedance_examples_obs <- function(df_station, station_name, threshold, colour = "#0072B2") {
  
  years <- select_example_years_obs(df_station, threshold)
  
  # Determine max y across valid years
  valid_years <- years[!is.na(unlist(years))]
  if(length(valid_years) == 0) {
    warning(paste("No valid example years for station:", station_name))
    return(NULL)
  }
  
  y_max <- df_station %>%
    filter(hydro_year %in% unlist(valid_years)) %>%
    summarise(max_val = max(rainfall_mm, na.rm = TRUE)) %>%
    pull(max_val) %>% {. + 10}  # add buffer
  
  # Create plots safely
  plots <- list()
  
  if(!is.na(years$muted)) {
    plots$muted <- plot_daily_example_obs(df_station, years$muted, threshold, paste0("Muted Year (", years$muted, ")"), y_max, colour)
  }
  if(!is.na(years$single)) {
    plots$single <- plot_daily_example_obs(df_station, years$single, threshold, paste0("Single Exceedance Year (", years$single, ")"), y_max, colour)
  }
  if(!is.na(years$high)) {
    plots$high <- plot_daily_example_obs(df_station, years$high, threshold, paste0("High Exceedance Year (", years$high, ")"), y_max, colour)
  }
  
  # Combine available plots
  wrap_plots(plots) +
    plot_annotation(title = station_name, theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
}

# Master function for example years per station
run_example_years_obs <- function(df_obs, station_name, output_dir) {
  
  if(!station_name %in% df_obs$station) {
    stop(paste("Station not found:", station_name))
  }
  
  df_station <- df_obs %>% filter(station == station_name)
  
  # Calculate thresholds
  thr_station <- calculate_rx1day_thresholds_obs(df_station)
  
  # Example-year plots (mid threshold)
  p_example_mid <- plot_exceedance_examples_obs(df_station, station_name, thr_station$thresholds["mid"])
  p_example_360 <- plot_exceedance_examples_obs(df_station, station_name, thr_station$thresholds["onein360"], colour = "#E69F00")
  
  # Save plots if they exist
  station_safe <- gsub(" ", "_", station_name)
  
  if(!is.null(p_example_mid)) {
    ggsave(
      filename = file.path(output_dir, paste0(station_safe, "_RX1day_ExampleYears_Mid.png")),
      plot = p_example_mid,
      width = fig_width_full,
      height = fig_height_tall,
      dpi = 300
    )
  }
  
  if(!is.null(p_example_360)) {
    ggsave(
      filename = file.path(output_dir, paste0(station_safe, "_RX1day_ExampleYears_360.png")),
      plot = p_example_360,
      width = fig_width_full,
      height = fig_height_tall,
      dpi = 300
    )
  }
  
  invisible(list(mid = p_example_mid, onein360 = p_example_360))
}
imap(region_output_dirs, function(region_dir, region_name) {
  
  # Get all stations in the region
  stations <- unique(combined_df_obs$station[combined_df_obs$region == region_name])
  
  # Loop through stations
  walk(stations, function(station_name) {
    
    # Subset data for this station
    df_station <- combined_df_obs %>% filter(station == station_name)
    
    # Calculate thresholds
    thr_station <- calculate_rx1day_thresholds_obs(df_station)
    
    # Check if station has any valid example years for either mid or 1-in-360 thresholds
    years_mid <- select_example_years_obs(df_station, thr_station$thresholds["mid"])
    years_360 <- select_example_years_obs(df_station, thr_station$thresholds["onein360"])
    
    all_years <- c(unlist(years_mid), unlist(years_360))
    
    # Skip station if no valid years
    if(all(is.na(all_years))) {
      message(paste("Skipping station (no valid example years):", station_name))
      return(NULL)
    }
    
    # Otherwise, run example-years analysis
    run_example_years_obs(combined_df_obs, station_name, region_dir)
  })
})


# boxplots ----------------------------------------------------------------

run_rx1day_boxplots_obs <- function(df_obs, station_name, output_dir) {
  
  if(!station_name %in% df_obs$station) {
    stop(paste("Station not found:", station_name))
  }
  
  df_station <- df_obs %>% filter(station == station_name)
  
  RX1day_df <- df_station %>%
    group_by(hydro_year) %>%
    summarise(RX1day = max(rainfall_mm, na.rm = TRUE), .groups = "drop") %>%
    rename(Year = hydro_year)
  
  thr_station <- calculate_rx1day_thresholds_obs(df_station)
  
  exc_mid <- count_exceedances_per_year_obs(df_station, thr_station$thresholds["mid"]) %>%
    rename(Year = hydro_year, exceedances = exceedance_days)
  
  exc_360 <- count_exceedances_per_year_obs(df_station, thr_station$thresholds["onein360"]) %>%
    rename(Year = hydro_year, exceedances = exceedance_days)
  
  df_mid  <- left_join(RX1day_df, exc_mid, by = "Year")
  df_360  <- left_join(RX1day_df, exc_360, by = "Year")
  
  rx1day_max <- max(RX1day_df$RX1day, na.rm = TRUE) * 1.1
  
  make_box <- function(df_plot, title) {
    
    exceed_max <- max(df_plot$exceedances, na.rm = TRUE)
    
    percent_per_bin <- df_plot %>%
      group_by(exceedances) %>%
      summarise(percent = n() / nrow(df_plot) * 100, .groups = "drop") %>%
      complete(exceedances = 0:exceed_max, fill = list(percent = 0))  # ensure all bins
    
    y_top <- max(df_plot$RX1day, na.rm = TRUE) * 1.12  # extra headroom
    
    ggplot(df_plot, aes(x = factor(exceedances), y = RX1day)) +
      
      geom_boxplot(
        fill = "grey85",
        colour = "black",
        alpha = 0.8,
        outlier.shape = NA
      ) +
      
      geom_jitter(
        data = function(d) {
          box_stats <- boxplot.stats(d$RX1day)$stats
          d %>% filter(RX1day > box_stats[5])
        },
        width = 0.15,
        alpha = 0.4,
        size = 1.2,
        colour = "black"
      ) +
      
      geom_text(
        data = percent_per_bin,
        aes(
          x = factor(exceedances),
          y = y_top,
          label = ifelse(percent < 0.1 & percent > 0, "<0.1%",
                         ifelse(percent == 0, "", paste0(round(percent, 1), "%")))
        ),
        inherit.aes = FALSE,
        size = 2.3,
        vjust = 0
      ) +
      
      theme_thesis +
      labs(
        x = "Number of exceedances per year",
        y = "RX1day (mm)",
        title = paste(station_name, "-", title)
      ) +
      guides(fill = "none") +
            scale_y_continuous(breaks = seq(0, y_top, by = 100)) +
      
      coord_cartesian(ylim = c(0, y_top), clip = "off") +
      scale_x_discrete(limits = as.character(0:exceed_max))
  }
  
  
  p_mid  <- make_box(df_mid,  "Mid RX1day Threshold")
  p_360  <- make_box(df_360, "1-in-360 Threshold")
  
  station_safe <- gsub(" ", "_", station_name)
  
  ggsave(file.path(output_dir, paste0(station_safe, "_RX1day_Boxplot_Mid.png")),
         p_mid, width = fig_width_full, height = fig_height_med, dpi = 300)
  
  ggsave(file.path(output_dir, paste0(station_safe, "_RX1day_Boxplot_360.png")),
         p_360, width = fig_width_full, height = fig_height_med, dpi = 300)
}

imap(region_output_dirs, function(region_dir, region_name) {
  
  stations <- unique(combined_df_obs$station[combined_df_obs$region == region_name])
  
  walk(stations, function(station_name) {
    
    run_RX1day_station_analysis_obs(combined_df_obs, station_name, region_dir)
    run_example_years_obs(combined_df_obs, station_name, region_dir)
    run_rx1day_boxplots_obs(combined_df_obs, station_name, region_dir)  
  })
})
