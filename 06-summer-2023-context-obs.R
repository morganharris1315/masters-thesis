# -------------------------------------------------------------------------
# 06-summer-2023-context-obs.R
# -------------------------------------------------------------------------
# Jan 2026
# Contextualise summer 2023 rainfall for observational stations.
# -------------------------------------------------------------------------

# Regions to process -------------------------------------------------------
regions_to_process <- c("coromandel", "far_north", "top_of_south", "waikato")

# Helper functions ---------------------------------------------------------
compute_RX1day_obs <- function(df_station, missing_day_threshold = 30) {
  df_station %>%
    group_by(hydro_year) %>%
    summarise(
      missing_days = sum(is.na(rainfall_mm)),
      RX1day_raw = if_else(all(is.na(rainfall_mm)), NA_real_, max(rainfall_mm, na.rm = TRUE)),
      RX1day = if_else(missing_days >= missing_day_threshold, NA_real_, RX1day_raw),
      .groups = "drop"
    ) %>%
    select(hydro_year, RX1day)
}

calculate_rx1day_thresholds_obs <- function(df_station_obs) {
  RX1day_df <- compute_RX1day_obs(df_station_obs)
  rx <- RX1day_df$RX1day

  # Keep single-threshold definition aligned with
  # 03a-blockmaxima-function-obs-singlethreshold.R:
  # threshold where 2/3 of annual RX1day values are above it.
  single_threshold <- quantile(
    rx,
    probs = 1 / 3,
    na.rm = TRUE,
    type = 7
  )

  list(
    thresholds = c(single = as.numeric(single_threshold))
  )
}

collapse_dates <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_character_)
  }
  paste(sort(unique(format(as.Date(x), "%Y-%m-%d"))), collapse = "; ")
}

# Create HY2023 + Summer 2023 example plots -------------------------------
create_example_year_plots <- function(df_station, station_name, threshold, output_dir, rx1day_hy2023_val = NA_real_) {
  hy_df <- df_station %>%
    filter(hydro_year == 2023) %>%
    arrange(observation_date)

  if (nrow(hy_df) == 0 || is.na(rx1day_hy2023_val)) {
    return(list(
      plot_generated = FALSE,
      plot_path = NA_character_
    ))
  }

  station_safe <- path_sanitize(station_name)
  plot_path <- file.path(output_dir, glue("{station_safe}_summer2023_context.png"))

  y_max <- max(hy_df$rainfall_mm, na.rm = TRUE)
  if (!is.finite(y_max)) y_max <- 10
  y_max <- y_max + 10

  p_hy <- ggplot(hy_df, aes(x = observation_date, y = rainfall_mm)) +
    geom_col(fill = "#0072B2", width= 1.5, na.rm = TRUE) +
    {if (is.finite(threshold)) geom_hline(yintercept = threshold, linetype = "dashed", colour = "#93acff", size = 1)} +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(
      title = glue("{station_name} — Hydrological Year 2023"),
      x = "Date",
      y = "Daily Rainfall (mm)"
    ) +
    theme_thesis

  ggsave(
    filename = plot_path,
    plot = p_hy,
    width = fig_width_standard,
    height = fig_height_standard,
    dpi = 300
  )

  list(
    plot_generated = TRUE,
    plot_path = plot_path
  )
}


# Region summary table -----------------------------------------------------
build_region_summer_2023_table <- function(df_obs, region_name, output_dir) {
  region_df <- df_obs %>%
    filter(region == region_name)

  stations <- sort(unique(region_df$station))

  summary_tbl <- map_dfr(stations, function(stn) {
    df_station <- region_df %>%
      filter(station == stn) %>%
      arrange(observation_date)

    thr <- calculate_rx1day_thresholds_obs(df_station)
    threshold_single <- thr$thresholds["single"]

    hy2023 <- df_station %>% filter(hydro_year == 2023)
    summer2023 <- hy2023 %>%
      filter(
        observation_date >= as.Date("2022-12-01"),
        observation_date <= as.Date("2023-02-28")
      )

    hy2023_total <- sum(hy2023$rainfall_mm, na.rm = TRUE)
    summer2023_total <- sum(summer2023$rainfall_mm, na.rm = TRUE)

    hy_exceed <- hy2023 %>% filter(rainfall_mm > threshold_single)
    summer_exceed <- summer2023 %>% filter(rainfall_mm > threshold_single)

    rx1day_all <- compute_RX1day_obs(df_station) %>% filter(!is.na(RX1day))
    rx1day_2023 <- hy2023 %>%
      filter(!is.na(rainfall_mm)) %>%
      slice_max(order_by = rainfall_mm, n = 1, with_ties = FALSE)

    rx1day_2023_val <- if (nrow(rx1day_2023) == 0) NA_real_ else rx1day_2023$rainfall_mm[[1]]
    rx1day_2023_date <- if (nrow(rx1day_2023) == 0) as.Date(NA) else rx1day_2023$observation_date[[1]]

    rx1day_percentile <- if (is.na(rx1day_2023_val) || nrow(rx1day_all) == 0) {
      NA_real_
    } else {
      round(mean(rx1day_all$RX1day <= rx1day_2023_val) * 100, 1)
    }

    plot_result <- create_example_year_plots(
      df_station = df_station,
      station_name = stn,
      threshold = threshold_single,
      output_dir = output_dir,
      rx1day_hy2023_val = rx1day_2023_val
    )

    tibble(
      region = region_name,
      station = stn,
      hydrological_year = 2023,
      annual_rainfall_mm = hy2023_total,
      summer_2023_rainfall_mm = summer2023_total,
      exceedance_threshold_mm = round(threshold_single, 2),
      exceedance_count_hy2023 = nrow(hy_exceed),
      exceedance_dates_hy2023 = collapse_dates(hy_exceed$observation_date),
      exceedance_count_summer2023 = nrow(summer_exceed),
      exceedance_dates_summer2023 = collapse_dates(summer_exceed$observation_date),
      rx1day_hy2023_mm = rx1day_2023_val,
      rx1day_hy2023_date = as.character(rx1day_2023_date),
      rx1day_hy2023_percentile = rx1day_percentile,
    )
  })

  write_csv(
    summary_tbl,
    file.path(output_dir, glue("{region_name}_summer_2023_station_summary.csv"))
  )

  summary_tbl
}

# Run for all selected regions --------------------------------------------
region_summaries <- map(
  regions_to_process,
  function(region_name) {
    region_context_dir <- glue("{base_raw_dir}/obs_data/{region_name}/summer_2023_context")
    build_region_summer_2023_table(
      df_obs = combined_df_obs,
      region_name = region_name,
      output_dir = region_context_dir
    )
  }
)

names(region_summaries) <- regions_to_process

summer_2023_all_regions_table <- bind_rows(region_summaries)

write_csv(
  summer_2023_all_regions_table,
  glue("{base_raw_dir}/obs_data/summer_2023_station_summary_all_regions.csv")
)

print(summer_2023_all_regions_table)


# Mapping -----------------------------------------------------------------


# Base Coromandel Plot ----------------------------------------------------

nz_map_df<-map_data('nz')

sites_path <- glue("{base_raw_dir}/Identified_Sites.csv")
identified_sites <- read.csv(sites_path)

coromandel_lat_log_data <-identified_sites %>% filter(CaseStudy == "Coromandel")
  
ggplot() +
  geom_polygon(
    data = nz_map_df,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    colour = "grey40",
    linewidth = 0.2
  ) +
  geom_point(
    data = coromandel_lat_log_data,
    aes(x = Longitude, y = Latitude),
    colour = "#0072B2",
    size = 2,
    alpha = 0.9
  ) +
  coord_quickmap(xlim = c(174.3, 176.2), ylim = c(-38.1, -35.9), expand = FALSE) +
  labs(
    title = "Coromandel stations",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_thesis



