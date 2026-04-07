# -------------------------------------------------------------------------
# 06-summer-2023-context-obs.R
# -------------------------------------------------------------------------
# Jan 2026
# Contextualise summer 2023 rainfall for observational stations.
# -------------------------------------------------------------------------

# Changes still to be made. 
# I want to add LiDAR data so that I have elevation on the map
# This should solve the issue of some of the dots going off the side of the land. 

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

format_day_month <- function(x) {
  x <- as.Date(x)
  if (is.na(x)) {
    return(NA_character_)
  }
  glue("{day(x)} {format(x, '%b')}")
}

format_day_month_year <- function(x) {
  x <- as.Date(x)
  if (is.na(x)) {
    return(NA_character_)
  }
  glue("{day(x)} {format(x, '%b %Y')}")
}

format_event_title <- function(start_date, end_date = start_date) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  if (is.na(start_date) || is.na(end_date)) {
    return(NA_character_)
  }
  
  if (start_date == end_date) {
    return(format_day_month_year(start_date))
  }
  
  if (year(start_date) == year(end_date) && month(start_date) == month(end_date)) {
    return(glue("{day(start_date)} to {day(end_date)} {format(end_date, '%b %Y')}"))
  }
  
  if (year(start_date) == year(end_date)) {
    return(glue("{format_day_month(start_date)} to {format_day_month_year(end_date)}"))
  }
  
  glue("{format_day_month_year(start_date)} to {format_day_month_year(end_date)}")
}

build_axis_breaks <- function(limits, n_breaks = 3) {
  if (length(limits) != 2 || any(!is.finite(limits))) {
    stop("limits must be a numeric vector of length 2 with finite values")
  }
  
  limits <- sort(limits)
  
  if (n_breaks <= 2 || diff(limits) == 0) {
    return(limits)
  }
  
  unique(round(seq(limits[1], limits[2], length.out = n_breaks), 1))
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
    geom_col(fill = "#93acff", width= 1.5, na.rm = TRUE) +
    {if (is.finite(threshold)) geom_hline(yintercept = threshold, linetype = "dashed", colour = "#93acff", size = 1)} +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(
      title = glue("{station_name} - 2023"),
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
  
  summary_tbl <- purrr::map_dfr(stations, function(stn) {
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
region_summaries <- purrr::map(
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

coromandel_lat_log_data <- identified_sites %>%
  filter(CaseStudy == "Coromandel")

p_coromandel_map <- ggplot() +
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
  coord_quickmap(xlim = c(175.2, 176.2), ylim = c(-37.5, -36.5), expand = FALSE) +
  scale_x_continuous(
    breaks = c(175.2, 175.7, 176.2),
    labels = scales::label_number(accuracy = 0.1),
    minor_breaks = NULL
  ) +
  labs(
    title = "Coromandel Stations"
  ) +
  theme_thesis

print(p_coromandel_map)

ggsave(file.path("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/obs_data/coromandel/coromandel_base_map.png"),
       plot = p_coromandel_map, width = 8, height = 7, dpi = 300)

coromandel_lat_log_data %>%
  summarise(
    n = n(),
    lon_min = min(Longitude, na.rm = TRUE),
    lon_max = max(Longitude, na.rm = TRUE),
    lat_min = min(Latitude, na.rm = TRUE),
    lat_max = max(Latitude, na.rm = TRUE)
  )


#  Coromandel key-event rainfall maps  ------------------------------

coromandel_clean_path <- glue("{base_raw_dir}/obs_data/cleaned_datasets/rain_coromandel_cleaned.csv")
coromandel_obs <- read_csv(coromandel_clean_path, show_col_types = FALSE) %>%
  mutate(
    observation_date = as.Date(observation_date),
    hydro_year = if_else(month(observation_date) >= 7,
                         year(observation_date) + 1,
                         year(observation_date))
  )

# Use the same station thresholds already used in the summary table so
# exceedance highlighting is directly consistent with reported values.
station_thresholds <- summer_2023_all_regions_table %>%
  filter(region == "coromandel") %>%
  transmute(
    station,
    threshold_single = as.numeric(exceedance_threshold_mm)
  ) %>%
  distinct(station, .keep_all = TRUE)

station_coords <- coromandel_obs %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  group_by(station) %>%
  summarise(
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop"
  )

event_definitions <- tribble(
  ~event_id, ~start_date, ~end_date,
  1L, as.Date("2022-07-25"), as.Date("2022-07-25"),
  2L, as.Date("2022-11-11"), as.Date("2022-11-11"),
  3L, as.Date("2022-12-14"), as.Date("2022-12-14"),
  4L, as.Date("2023-01-09"), as.Date("2023-01-11"),
  5L, as.Date("2023-01-27"), as.Date("2023-01-28"),
  6L, as.Date("2023-02-12"), as.Date("2023-02-14")
) %>%
  mutate(event_title = map2_chr(start_date, end_date, format_event_title))

chiltern_site <- tibble(
  station = "Chiltern",
  latitude = -36.81804,
  longitude = 175.53654
)

build_event_map <- function(event_row, base_map_df, xlim = c(175.2, 176.2), ylim = c(-37.6, -36.4)) {
  legend_labels <- c(
    "Daily Rainfall (mm)",
    "Above Heavy Threshold"
  )
  
  event_data <- coromandel_obs %>%
    filter(
      observation_date >= event_row$start_date,
      observation_date <= event_row$end_date
    ) %>%
    group_by(station) %>%
    filter(!is.na(rainfall_mm)) %>%
    slice_max(order_by = rainfall_mm, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      station,
      event_rainfall_mm = rainfall_mm,
      event_date_used = observation_date
    ) %>%
    left_join(station_coords, by = "station") %>%
    left_join(station_thresholds, by = "station") %>%
    mutate(is_above_threshold = !is.na(threshold_single) & event_rainfall_mm > threshold_single) %>%
    filter(!is.na(latitude), !is.na(longitude), !is.na(event_rainfall_mm))
  
  message(
    glue(
      "Event {event_row$event_title}: {sum(event_data$is_above_threshold, na.rm = TRUE)} stations above threshold"
    )
  )
  
  ggplot() +
    geom_polygon(
      data = base_map_df,
      aes(x = long, y = lat, group = group),
      fill = "grey92",
      colour = "grey40",
      linewidth = 0.2
    ) +
    geom_point(
      data = event_data,
      aes(
        x = longitude,
        y = latitude,
        fill = legend_labels[1],
        colour = legend_labels[1]
      ),
      shape = 21,
      size = 2,
      stroke = 0.4,
      alpha = 0.95
    ) +
    geom_point(
      data = event_data %>% filter(is_above_threshold),
      aes(
        x = longitude,
        y = latitude,
        fill = legend_labels[2],
        colour = legend_labels[2]
      ),
      shape = 21,
      size = 3,
      stroke = 1.1
    ) +
    geom_text(
      data = chiltern_site,
      aes(x = longitude, y = latitude, label = station),
      nudge_x = -0.05,
      hjust = 1,
      vjust = 0.5,
      size = 2.6,
      colour = "red",
      fontface = "bold"
    ) +
    geom_text(
      data = event_data,
      aes(x = longitude, y = latitude, label = round(event_rainfall_mm, 0)),
      nudge_y = 0.045,
      size = 2.5,
      colour = "grey10",
      check_overlap = TRUE
    ) +
    coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_x_continuous(
      breaks = c(xlim[1], mean(xlim), xlim[2]),
      labels = scales::label_number(accuracy = 0.1),
      minor_breaks = NULL
    ) +
    scale_y_continuous(
      breaks = seq(ylim[1], ylim[2], by = 0.4),
      minor_breaks = NULL
    ) +
    labs(
      title = event_row$event_title,
      x = NULL,
      y = NULL
    ) +
    scale_colour_manual(
      values = c(
        "Daily Rainfall (mm)" = "grey20",
        "Above Heavy Threshold" = "#93acff"
      ),
      breaks = legend_labels,
      name = NULL
    ) +
    scale_fill_manual(
      values = c(
        "Daily Rainfall (mm)" = "grey45",
        "Above Heavy Threshold" = "grey45"
      ),
      breaks = legend_labels,
      name = NULL
    ) +
    guides(
      fill = "none",
      colour = guide_legend(
        override.aes = list(
          shape = 21,
          size = c(2.7, 3.4),
          stroke = c(0.4, 1.1),
          fill = c("grey45", "grey45"),
          colour = c("grey20", "#93acff"),
          alpha = 1
        )
      )
    ) +
    theme_thesis +
    theme(
      axis.title = element_blank(),
      legend.position = "none"
    )
}

event_maps <- purrr::map(
  split(event_definitions, event_definitions$event_id),
  ~ build_event_map(.x, nz_map_df)
)

event_map_grid <- patchwork::wrap_plots(event_maps, ncol = 3, nrow = 2)

# Draw a single shared legend for the full panel on the right-hand side.
p_coromandel_event_panel <- event_map_grid +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(p_coromandel_event_panel)

ggsave(
  filename = glue("{base_raw_dir}/obs_data/coromandel/coromandel_key_event_maps_3x2.png"),
  plot = p_coromandel_event_panel,
  dpi = 300)


# Adding Lidar Data -------------------------------------------------------
#install.packages("lidR")
#install.packages("terra")
#install.packages("rayshader")
#install.packages("terrainr")

# Load the libraries
library(lidR)
library(terra)
library(rayshader)
library(viridis)
library(terrainr)


# Build Coromandel LiDAR hillshade basemap ---------------------------------
lidar_dir <- glue("{base_raw_dir}/obs_data/coromandel/LiDAR data")
lidar_tile_paths <- dir_ls(lidar_dir, regexp = "\\.tif$", recurse = FALSE)

if (length(lidar_tile_paths) == 0) {
  stop(glue("No .tif files found in LiDAR directory: {lidar_dir}"))
}

lidar_tiles <- lapply(lidar_tile_paths, terra::rast)
lidar_mosaic <- if (length(lidar_tiles) == 1) {
  lidar_tiles[[1]]
} else {
  terra::mosaic(terra::sprc(lidar_tiles))
}

# Ensure station lon/lat and NZ outline can be plotted on the same CRS.
lidar_wgs84 <- terra::project(lidar_mosaic, "EPSG:4326")

coromandel_extent <- terra::ext(175.2, 176.2, -37.5, -36.5)
lidar_wgs84_crop <- terra::crop(lidar_wgs84, coromandel_extent)

# Build simple grayscale hillshade background.
lidar_slope <- terra::terrain(lidar_wgs84_crop, v = "slope", unit = "radians")
lidar_aspect <- terra::terrain(lidar_wgs84_crop, v = "aspect", unit = "radians")
lidar_hillshade <- terra::shade(lidar_slope, lidar_aspect, angle = 45, direction = 315)
names(lidar_hillshade) <- "hillshade"

lidar_hillshade_df <- as.data.frame(lidar_hillshade, xy = TRUE, na.rm = TRUE)

p_coromandel_lidar_base_map <- ggplot() +
  geom_raster(
    data = lidar_hillshade_df,
    aes(x = x, y = y, fill = hillshade)
  ) +
  scale_fill_gradient(
    low = "grey15",
    high = "grey92",
    name = "Hillshade"
  ) +
  geom_polygon(
    data = nz_map_df,
    aes(x = long, y = lat, group = group),
    fill = NA,
    colour = "grey20",
    linewidth = 0.35
  ) +
  geom_point(
    data = coromandel_lat_log_data,
    aes(x = Longitude, y = Latitude),
    shape = 21,
    fill = "grey55",
    colour = "black",
    size = 2
  ) +
  coord_quickmap(
    xlim = c(175.2, 176.2),
    ylim = c(-37.5, -36.5),
    expand = FALSE
  ) +
  labs(
    title = "Coromandel Stations with LiDAR Hillshade Background",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_thesis +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "right"
  )

print(p_coromandel_lidar_base_map)

ggsave(
  filename = glue("{base_raw_dir}/obs_data/coromandel/coromandel_base_map_lidar_hillshade.png"),
  plot = p_coromandel_lidar_base_map,
  width = 8,
  height = 7,
  dpi = 300
)
