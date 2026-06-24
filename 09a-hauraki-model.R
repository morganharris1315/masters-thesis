# -------------------------------------------------------------------------
# 09a-hauraki-model.R
# -------------------------------------------------------------------------
# Hauraki version of figure 2.
# -------------------------------------------------------------------------

# Packages ----------------------------------------------------------------
library(RNetCDF)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(purrr)
library(readr)
library(lubridate)
library(maps)
library(scales)
library(terra)
library(sf)

# Style constants ----------------------------------------------------------
heavy_col <- "lightslateblue"
extreme_col <- "darkslateblue"
box_colour_light <- "lavender"
box_colour_dark <- "lightslateblue"

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

hauraki_current_dir <- file.path(weatherathome_dir, "current_decade")
hauraki_future_dir  <- file.path(weatherathome_dir, "3k_warmer")

hauraki_figure1_file <- file.path(base_raw_dir, "obs_data", "Figure1_hauraki_context.png")
hauraki_figure2_file <- file.path(weatherathome_dir, "Figure2_hauraki_cell_2x3.png")
hauraki_figure5_file <- file.path(weatherathome_dir, "Figure5_hauraki_ge5_heatmap.png")

# Inputs -------------------------------------------------------------------
hauraki_matched_cell <- data.frame(lon_index = 30L, lat_index = 16L)

# Shared helpers -------------------------------------
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

pct_change <- function(future_value, current_value) {
  if (isTRUE(all.equal(current_value, 0))) return(NA_real_)
  100 * (future_value - current_value) / current_value
}


# Figure 2 Hauraki --------------------------------------------------------

hauraki_current_files <- list_nc_files(hauraki_current_dir)
hauraki_future_files  <- list_nc_files(hauraki_future_dir)

cd_df_hauraki <- build_annual_cell_df(
  hauraki_current_files,
  hauraki_matched_cell$lon_index,
  hauraki_matched_cell$lat_index,
  "Current Day"
)

fp_df_hauraki <- build_annual_cell_df(
  hauraki_future_files,
  hauraki_matched_cell$lon_index,
  hauraki_matched_cell$lat_index,
  "Future Projection"
)

heavy_cd_hauraki <- as.numeric(quantile(cd_df_hauraki$RX1day, 0.33, na.rm = TRUE, type= 7))
heavy_fp_hauraki <- as.numeric(quantile(fp_df_hauraki$RX1day, 0.33, na.rm = TRUE, type= 7))
extreme_cd_hauraki <- as.numeric(quantile(cd_df_hauraki$RX1day, 0.90, na.rm = TRUE, type= 7))
extreme_fp_hauraki <- as.numeric(quantile(fp_df_hauraki$RX1day, 0.90, na.rm = TRUE, type= 7))

cd_df_hauraki <- cd_df_hauraki %>%
  mutate(
    heavy_threshold = heavy_cd_hauraki,
    extreme_threshold = extreme_cd_hauraki,
    heavy_days = count_heavy_days(
      hauraki_current_files,
      hauraki_matched_cell$lon_index,
      hauraki_matched_cell$lat_index,
      heavy_cd_hauraki
    )
  )

fp_df_hauraki <- fp_df_hauraki %>%
  mutate(
    heavy_threshold = heavy_fp_hauraki,
    extreme_threshold = extreme_fp_hauraki,
    heavy_days = count_heavy_days(
      hauraki_future_files,
      hauraki_matched_cell$lon_index,
      hauraki_matched_cell$lat_index,
      heavy_fp_hauraki
    )
  )

fp_df_cd_thresholds_hauraki <- fp_df_hauraki %>%
  mutate(
    heavy_threshold = heavy_cd_hauraki,
    extreme_threshold = extreme_cd_hauraki,
    heavy_days = count_heavy_days(
      hauraki_future_files,
      hauraki_matched_cell$lon_index,
      hauraki_matched_cell$lat_index,
      heavy_cd_hauraki
    )
  )

# Derived outputs 
hist_cd_hauraki <- build_histogram_df(cd_df_hauraki$heavy_days)
hist_fp_hauraki <- build_histogram_df(fp_df_cd_thresholds_hauraki$heavy_days)

heavy_days_max_hauraki <- max(
  c(cd_df_hauraki$heavy_days, fp_df_cd_thresholds_hauraki$heavy_days),
  na.rm = TRUE
)

heavy_change_pct_hauraki <- pct_change(heavy_fp_hauraki, heavy_cd_hauraki)
extreme_change_pct_hauraki <- pct_change(extreme_fp_hauraki, extreme_cd_hauraki)

rx1day_upper_shared_hauraki <- max(
  c(cd_df_hauraki$RX1day, fp_df_hauraki$RX1day),
  na.rm = TRUE
) * 1.05

# Plots
p2a_hauraki <- plot_rx1day_ts(cd_df_hauraki, "(a)", FALSE, y_upper = rx1day_upper_shared_hauraki)
p2b_hauraki <- plot_rx1day_ts(
  fp_df_hauraki,
  "(b)",
  TRUE,
  heavy_change = heavy_change_pct_hauraki,
  extreme_change = extreme_change_pct_hauraki,
  y_upper = rx1day_upper_shared_hauraki
)

p2c_hauraki <- plot_heavy_hist(hist_cd_hauraki, "(c)", heavy_days_max_hauraki)
p2d_hauraki <- plot_heavy_hist(hist_fp_hauraki, "(d)", heavy_days_max_hauraki)
p2e_hauraki <- plot_quadrant_heatmap(cd_df_hauraki, "(e)", 4L, heavy_days_max_hauraki)
p2f_hauraki <- plot_quadrant_heatmap(fp_df_cd_thresholds_hauraki, "(f)", 4L, heavy_days_max_hauraki)

col_left_hauraki <- make_column_header("Current Day") /
  p2a_hauraki /
  p2c_hauraki /
  p2e_hauraki +
  plot_layout(heights = c(0.09, 1, 1, 1))

col_right_hauraki <- make_column_header("Future Projection") /
  p2b_hauraki /
  p2d_hauraki /
  p2f_hauraki +
  plot_layout(heights = c(0.09, 1, 1, 1))

figure2_hauraki_plot <- col_left_hauraki | col_right_hauraki

# Figure 5 (Hauraki)
p2e5_hauraki <- plot_quadrant_heatmap(cd_df_hauraki, "(a)", 5L, heavy_days_max_hauraki)
p2f5_hauraki <- plot_quadrant_heatmap(fp_df_cd_thresholds_hauraki, "(b)", 5L, heavy_days_max_hauraki)

figure5_hauraki_plot <- (make_column_header("Current Day") / p2e5_hauraki) |
  (make_column_header("Future Projection") / p2f5_hauraki)

# Save --------------------------------------------------------------------
ggsave(hauraki_figure2_file, figure2_hauraki_plot, width = 7, height = 8, dpi = 2000)
ggsave(hauraki_figure5_file, figure5_hauraki_plot, width = 7, height = 2.6, dpi = 2000)
