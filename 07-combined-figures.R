# -------------------------------------------------------------------------
# 07-combined-figures.R
# -------------------------------------------------------------------------
# Mar 2026
# Build combined paper figures:
# - Figure 1: observational context + Chiltern station panels
# - Figure 2: matched weather@home grid-cell six-panel workflow
# -------------------------------------------------------------------------

# Packages ----------------------------------------------------------------
library(RNetCDF)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(lubridate)

# Style constants ----------------------------------------------------------
heavy_col <- "#93acff"
extreme_col <- "#5f66DB"
box_colour_light <- "#eff2ff"
box_colour_dark <- "#6f8dff"

if (!exists("theme_thesis")) {
  theme_thesis <- theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7),
      strip.text = element_text(size = 8, face = "bold"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2),
      panel.spacing = unit(0.5, "lines")
    )
}

theme_model_axes <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = "black", linewidth = 0.3),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)
)

# Paths --------------------------------------------------------------------
thesis_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis"
base_raw_dir <- file.path(thesis_dir, "Compound Events")
weatherathome_dir <- file.path(base_raw_dir, "model_data")

current_dir <- file.path(weatherathome_dir, "current_decade")
future_dir <- file.path(weatherathome_dir, "3k_warmer")
chiltern_obs_file <- file.path(base_raw_dir, "obs_data", "cleaned_datasets", "rain_coromandel_cleaned.csv")

figure1_file <- file.path(base_raw_dir, "obs_data", "Figure1_combined.png")
figure2_file <- file.path(weatherathome_dir, "Figure2_coromandel_cell_2x3.png")

# Inputs -------------------------------------------------------------------
matched_cell <- data.frame(lon_index = 30L, lat_index = 16L)
chiltern_site <- tibble(
  station = "Chiltern",
  latitude = -36.81804,
  longitude = 175.53654
)

# Shared helpers -----------------------------------------------------------
list_nc_files <- function(directory_path) {
  files <- list.files(directory_path, pattern = "\\.nc$", full.names = TRUE)
  if (length(files) == 0) stop("No .nc files found in: ", directory_path)
  files
}

extract_cell_daily_series <- function(file_path, lon_index, lat_index, var_name = "item5216_daily_mean") {
  nc <- open.nc(file_path)
  on.exit(close.nc(nc), add = TRUE)
  
  pr <- var.get.nc(nc, var_name)
  pr_dim <- dim(pr)
  
  if (length(pr_dim) == 4 && pr_dim[3] == 1) {
    pr <- pr[, , 1, ]
    pr_dim <- dim(pr)
  }
  
  if (length(pr_dim) != 3) stop("Unexpected precipitation dimensions in file: ", basename(file_path))
  
  if (lon_index < 1 || lon_index > pr_dim[1] || lat_index < 1 || lat_index > pr_dim[2]) {
    stop(
      "Requested lon/lat index is outside file grid range. ",
      "lon_index=", lon_index, " (max ", pr_dim[1], "), ",
      "lat_index=", lat_index, " (max ", pr_dim[2], ")."
    )
  }
  
  as.numeric(pr[lon_index, lat_index, ])
}

build_annual_cell_df <- function(nc_files, lon_index, lat_index, period_label) {
  bind_rows(lapply(seq_along(nc_files), function(i) {
    daily <- extract_cell_daily_series(nc_files[i], lon_index, lat_index)
    data.frame(
      Year = i,
      RX1day = max(daily, na.rm = TRUE),
      file = basename(nc_files[i]),
      Period = period_label,
      stringsAsFactors = FALSE
    )
  }))
}

count_heavy_days <- function(nc_files, lon_index, lat_index, heavy_threshold) {
  vapply(nc_files, function(f) {
    daily <- extract_cell_daily_series(f, lon_index, lat_index)
    sum(daily > heavy_threshold, na.rm = TRUE)
  }, numeric(1))
}

calc_cumulative_props <- function(exceed_vec) {
  max_k <- max(exceed_vec, na.rm = TRUE)
  n_years <- sum(!is.na(exceed_vec))
  
  data.frame(
    days = 0:max_k,
    cum_prop_years = sapply(0:max_k, function(k) sum(exceed_vec >= k, na.rm = TRUE) / n_years)
  )
}

build_histogram_df <- function(heavy_days_vec) {
  cumulative_df <- calc_cumulative_props(heavy_days_vec)
  
  data.frame(days = heavy_days_vec) %>%
    count(days, name = "n") %>%
    complete(days = 0:max(days, na.rm = TRUE), fill = list(n = 0)) %>%
    mutate(prop_years = n / sum(n)) %>%
    left_join(cumulative_df, by = "days") %>%
    mutate(
      cumulative_label = case_when(
        is.na(cum_prop_years) ~ NA_character_,
        cum_prop_years > 0 & (cum_prop_years * 100) < 1 ~ "<1%",
        TRUE ~ paste0(round(cum_prop_years * 100, 0), "%")
      )
    )
}

pct_change <- function(future_value, current_value) {
  if (isTRUE(all.equal(current_value, 0))) return(NA_real_)
  100 * (future_value - current_value) / current_value
}

plot_rx1day_ts <- function(df_period, panel_tag, include_change = FALSE, heavy_change = NA_real_, extreme_change = NA_real_, show_points = TRUE) {
  heavy_label <- if (!include_change) {
    sprintf("Heavy %.1f mm", unique(df_period$heavy_threshold))
  } else {
    sprintf("Heavy %.1f mm, %+0.1f%%", unique(df_period$heavy_threshold), heavy_change)
  }
  
  extreme_label <- if (!include_change) {
    sprintf("Extreme %.1f mm", unique(df_period$extreme_threshold))
  } else {
    sprintf("Extreme %.1f mm, %+0.1f%%", unique(df_period$extreme_threshold), extreme_change)
  }
  
  y_max <- max(df_period$RX1day, na.rm = TRUE) * 1.05
  
  p <- ggplot(df_period, aes(x = Year, y = RX1day)) +
    geom_line(colour = "black", linewidth = 0.35) +
    geom_hline(aes(yintercept = heavy_threshold), colour = heavy_col, linetype = "solid", linewidth = 0.9) +
    geom_hline(aes(yintercept = extreme_threshold), colour = extreme_col, linetype = "solid", linewidth = 0.9) +
    annotate("text", x = Inf, y = y_max, label = heavy_label, hjust = 1.03, vjust = 1.6, size = 2.7, colour = heavy_col) +
    annotate("text", x = Inf, y = y_max, label = extreme_label, hjust = 1.03, vjust = 3.1, size = 2.7, colour = extreme_col) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(title = panel_tag, x = "Year", y = "RX1day (mm)") +
    theme_thesis +
    theme_model_axes
  
  if (isTRUE(show_points)) {
    p <- p + geom_point(colour = "black", size = 0.45)
  }
  
  p
}

plot_heavy_hist <- function(hist_df, panel_tag) {
  max_exceed <- max(hist_df$days, na.rm = TRUE)
  
  ggplot(hist_df, aes(x = days, y = prop_years)) +
    geom_col(width = 0.9, fill = heavy_col) +
    geom_text(aes(label = cumulative_label), vjust = -0.35, size = 2.3, na.rm = TRUE) +
    scale_x_continuous(breaks = 0:max_exceed, limits = c(-0.5, max_exceed + 0.5), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = panel_tag, x = "Number of heavy days", y = "Proportion of years") +
    coord_cartesian(clip = "off") +
    theme_thesis +
    theme_model_axes +
    theme(
      panel.grid.major = element_line(colour = "grey82", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
      axis.line = element_line(colour = "black", linewidth = 0.2),
      plot.margin = margin(8, 10, 18, 8)
    )
}

build_quadrant_df <- function(df_period, heavy_cutoff = 4L) {
  exceed_max <- max(df_period$heavy_days, na.rm = TRUE)
  
  tile_df <- df_period %>%
    mutate(
      exceed_group = if_else(heavy_days < heavy_cutoff, paste0("<", heavy_cutoff, " heavy days"), paste0(">=", heavy_cutoff, " heavy days")),
      rx_group = if_else(RX1day < extreme_threshold, "below Extreme", "at/above Extreme")
    ) %>%
    count(exceed_group, rx_group, name = "n_years") %>%
    complete(
      exceed_group = c(paste0("<", heavy_cutoff, " heavy days"), paste0(">=", heavy_cutoff, " heavy days")),
      rx_group = c("below Extreme", "at/above Extreme"),
      fill = list(n_years = 0)
    ) %>%
    mutate(
      pct_years = 100 * n_years / sum(n_years),
      pct_label = paste0(round(pct_years, 1), "%"),
      xmin = if_else(exceed_group == paste0("<", heavy_cutoff, " heavy days"), -0.5, heavy_cutoff),
      xmax = if_else(exceed_group == paste0("<", heavy_cutoff, " heavy days"), heavy_cutoff, exceed_max + 0.5),
      ymin = if_else(rx_group == "below Extreme", 0, unique(df_period$extreme_threshold)),
      ymax = if_else(rx_group == "below Extreme", unique(df_period$extreme_threshold), max(df_period$RX1day, na.rm = TRUE) * 1.02),
      xmid = (xmin + xmax) / 2,
      ymid = (ymin + ymax) / 2
    )
  
  list(tile_df = tile_df, exceed_max = exceed_max)
}

plot_quadrant_heatmap <- function(df_period, panel_tag, heavy_cutoff = 4L) {
  hm <- build_quadrant_df(df_period, heavy_cutoff = heavy_cutoff)
  
  ggplot(hm$tile_df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = pct_years), colour = "white", linewidth = 0.3) +
    geom_text(aes(x = xmid, y = ymid, label = pct_label), colour = "black", size = 3.0, fontface = "bold") +
    geom_vline(xintercept = heavy_cutoff, colour = "black", linewidth = 0.35) +
    geom_hline(yintercept = unique(df_period$extreme_threshold), colour = "black", linewidth = 0.35) +
    scale_x_continuous(breaks = 0:hm$exceed_max, limits = c(-0.5, hm$exceed_max + 0.5), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, max(df_period$RX1day, na.rm = TRUE) * 1.02), expand = expansion(mult = c(0, 0))) +
    scale_fill_gradient(low = box_colour_light, high = box_colour_dark, limits = c(0, 100), guide = "none") +
    labs(title = panel_tag, x = "Number of heavy days", y = "RX1day (mm)") +
    theme_thesis +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.35),
      axis.line = element_blank()
    )
}

make_column_header <- function(label_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label_text, fontface = "bold", size = 3.5) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

# Figure 1 -----------------------------------------------------------------
if (!exists("p_coromandel_event_panel")) {
  stop("Object 'p_coromandel_event_panel' not found. Run 06-summer-2023-context-obs.R first in this R session.")
}

if (!file.exists(chiltern_obs_file)) {
  stop("Chiltern observational dataset not found: ", chiltern_obs_file)
}

chiltern_obs <- readr::read_csv(chiltern_obs_file, show_col_types = FALSE) %>%
  mutate(
    observation_date = as.Date(observation_date),
    hydro_year = if_else(month(observation_date) >= 7, year(observation_date) + 1L, year(observation_date))
  ) %>%
  filter(station == chiltern_site$station) %>%
  arrange(observation_date)

if (nrow(chiltern_obs) == 0) {
  stop("No rows found for station '", chiltern_site$station, "' in ", chiltern_obs_file)
}

chiltern_rx <- chiltern_obs %>%
  group_by(hydro_year) %>%
  summarise(RX1day = max(rainfall_mm, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(RX1day))

heavy_obs <- as.numeric(quantile(chiltern_rx$RX1day, probs = 0.33, na.rm = TRUE, type = 7))
extreme_obs <- as.numeric(quantile(chiltern_rx$RX1day, probs = 0.90, na.rm = TRUE, type = 7))

chiltern_rx_plot <- chiltern_rx %>%
  mutate(
    Year = hydro_year,
    heavy_threshold = heavy_obs,
    extreme_threshold = extreme_obs
  )

p1b <- plot_rx1day_ts(chiltern_rx_plot, panel_tag = "(b)", include_change = FALSE, show_points = TRUE) +
  geom_point(data = subset(chiltern_rx_plot, hydro_year == 2023), aes(x = Year, y = RX1day), colour = "red3", size = 1.7) +
  geom_text(
    data = subset(chiltern_rx_plot, hydro_year == 2023),
    aes(x = Year, y = RX1day, label = "2023"),
    nudge_y = max(chiltern_rx_plot$RX1day, na.rm = TRUE) * 0.03,
    colour = "red3",
    size = 2.8,
    fontface = "bold"
  )

hy2023_df <- chiltern_obs %>%
  filter(hydro_year == 2023) %>%
  arrange(observation_date)

if (nrow(hy2023_df) == 0) {
  stop("No data found for Chiltern hydrological year 2023.")
}

y_max_2023 <- max(hy2023_df$rainfall_mm, na.rm = TRUE)
if (!is.finite(y_max_2023)) y_max_2023 <- 10

axis_breaks <- seq(min(hy2023_df$observation_date, na.rm = TRUE), max(hy2023_df$observation_date, na.rm = TRUE), by = "1 month")

p1c <- ggplot(hy2023_df, aes(x = observation_date, y = rainfall_mm)) +
  geom_col(fill = "black", width = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = heavy_obs, linetype = "solid", colour = heavy_col, linewidth = 1) +
  annotate("text", x = Inf, y = y_max_2023 * 1.02, label = sprintf("Heavy %.1f mm", heavy_obs), hjust = 1.03, vjust = 1.5, size = 2.7, colour = heavy_col) +
  scale_x_date(breaks = axis_breaks, date_labels = "%b") +
  scale_y_continuous(limits = c(0, 300)) +
  labs(title = "(c)", x = "Date", y = "Daily Rainfall (mm)") +
  theme_thesis +
  theme_model_axes

# Add (a) label only, no title change to event panel internals
p1a <- p_coromandel_event_panel +
  plot_annotation(title = "(a)") &
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5))

figure1_plot <- p1a / (p1b | p1c) +
  plot_layout(heights = c(1, 1))

figure1_plot

# Figure 2 -----------------------------------------------------------------
current_files <- list_nc_files(current_dir)
future_files <- list_nc_files(future_dir)

cd_df <- build_annual_cell_df(current_files, matched_cell$lon_index, matched_cell$lat_index, "Current Day")
fp_df <- build_annual_cell_df(future_files, matched_cell$lon_index, matched_cell$lat_index, "Future Projection")

heavy_cd <- as.numeric(quantile(cd_df$RX1day, probs = 0.33, na.rm = TRUE, type = 7))
heavy_fp <- as.numeric(quantile(fp_df$RX1day, probs = 0.33, na.rm = TRUE, type = 7))
extreme_cd <- as.numeric(quantile(cd_df$RX1day, probs = 0.90, na.rm = TRUE, type = 7))
extreme_fp <- as.numeric(quantile(fp_df$RX1day, probs = 0.90, na.rm = TRUE, type = 7))

cd_df <- cd_df %>%
  mutate(
    heavy_threshold = heavy_cd,
    extreme_threshold = extreme_cd,
    heavy_days = count_heavy_days(current_files, matched_cell$lon_index, matched_cell$lat_index, heavy_cd)
  )

fp_df <- fp_df %>%
  mutate(
    heavy_threshold = heavy_fp,
    extreme_threshold = extreme_fp,
    heavy_days = count_heavy_days(future_files, matched_cell$lon_index, matched_cell$lat_index, heavy_fp)
  )

hist_cd <- build_histogram_df(cd_df$heavy_days)
hist_fp <- build_histogram_df(fp_df$heavy_days)

heavy_change_pct <- pct_change(heavy_fp, heavy_cd)
extreme_change_pct <- pct_change(extreme_fp, extreme_cd)

p2a <- plot_rx1day_ts(cd_df, panel_tag = "(a)", include_change = FALSE)
p2b <- plot_rx1day_ts(fp_df, panel_tag = "(b)", include_change = TRUE, heavy_change = heavy_change_pct, extreme_change = extreme_change_pct)
p2c <- plot_heavy_hist(hist_cd, panel_tag = "(c)")
p2d <- plot_heavy_hist(hist_fp, panel_tag = "(d)")
p2e <- plot_quadrant_heatmap(cd_df, panel_tag = "(e)", heavy_cutoff = 4L)
p2f <- plot_quadrant_heatmap(fp_df, panel_tag = "(f)", heavy_cutoff = 4L)

col_left <- make_column_header("Current Day") / p2a / p2c / p2e + plot_layout(heights = c(0.09, 1, 1, 1))
col_right <- make_column_header("Future Projection") / p2b / p2d / p2f + plot_layout(heights = c(0.09, 1, 1, 1))

figure2_plot <- col_left | col_right

# Save outputs -------------------------------------------------------------
ggsave(filename = figure1_file, plot = figure1_plot, width = 11, height = 8.6, dpi = 300)
ggsave(filename = figure2_file, plot = figure2_plot, width = 11, height = 10.4, dpi = 300)

message("Figure 1 saved to: ", figure1_file)
message("Figure 2 saved to: ", figure2_file)
message("Matched cell used for Figure 2: lon_index=", matched_cell$lon_index, ", lat_index=", matched_cell$lat_index)
message(sprintf("Figure 2 current thresholds -> Heavy: %.2f mm, Extreme: %.2f mm", heavy_cd, extreme_cd))
message(sprintf("Figure 2 future thresholds  -> Heavy: %.2f mm (%+.1f%%), Extreme: %.2f mm (%+.1f%%)", heavy_fp, heavy_change_pct, extreme_fp, extreme_change_pct))
