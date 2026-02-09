# -------------------------------------------------------------------------
# 03b-blockmaxima-function-modeldata.R
# -------------------------------------------------------------------------
# Jan 2026
# RX1day processing for model data using a single threshold definition
# Threshold: value where 2/3 of annual RX1day values are above it
# -------------------------------------------------------------------------

# Read raw data ------------------------------------------------------------
waikato_CD_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/WaikatoCurrentClimate.csv")
waikato_FP_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/Waikato3DegClimate.csv")
napier_CD_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/NapierCurrentClimate.csv")
napier_FP_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/Napier3DegClimate.csv")
northland_CD_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/NorthlandCurrentClimate.csv")
northland_FP_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/Northland3DegClimate.csv")
milford_CD_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/wah_ens_daily_rain_MilfordSound_CD.csv")
milford_FP_raw <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/wah_ens_daily_rain_MilfordSound_3deg.csv")

# Prepare model data -------------------------------------------------------
prepare_model_data <- function(df) {
  colnames(df) <- as.character(1:360)
  df %>%
    mutate(
      RX1day = apply(across(1:360), 1, max, na.rm = TRUE),
      Year = row_number()
    )
}

waikato_CD <- prepare_model_data(waikato_CD_raw)
waikato_FP <- prepare_model_data(waikato_FP_raw)
napier_CD <- prepare_model_data(napier_CD_raw)
napier_FP <- prepare_model_data(napier_FP_raw)
northland_CD <- prepare_model_data(northland_CD_raw)
northland_FP <- prepare_model_data(northland_FP_raw)
milford_CD <- prepare_model_data(milford_CD_raw)
milford_FP <- prepare_model_data(milford_FP_raw)

regions_mod <- c("waikato", "napier", "northland", "milford")
region_labels <- c(
  waikato = "Waikato",
  napier = "Napier",
  northland = "Northland",
  milford = "Milford Sound"
)
period_labels <- c(CD = "Current Day", FP = "Future Projection")
box_colour <- "#4F46E5"

# Calculate thresholds -----------------------------------------------------
calculate_rx1day_threshold <- function(df) {
  rx <- df$RX1day
  threshold <- as.numeric(quantile(rx, probs = 1 / 3, na.rm = TRUE, type = 7))

  daily_vec <- df %>%
    select(1:360) %>%
    as.matrix() %>%
    as.numeric()

  n_days <- sum(!is.na(daily_vec))
  prop_exceed <- sum(daily_vec > threshold, na.rm = TRUE) / n_days

  list(
    threshold = threshold,
    proportion = prop_exceed
  )
}

thr_list <- list(
  waikato_CD = calculate_rx1day_threshold(waikato_CD),
  waikato_FP = calculate_rx1day_threshold(waikato_FP),
  napier_CD = calculate_rx1day_threshold(napier_CD),
  napier_FP = calculate_rx1day_threshold(napier_FP),
  northland_CD = calculate_rx1day_threshold(northland_CD),
  northland_FP = calculate_rx1day_threshold(northland_FP),
  milford_CD = calculate_rx1day_threshold(milford_CD),
  milford_FP = calculate_rx1day_threshold(milford_FP)
)

threshold_label <- function(thr) {
  paste0(
    "2/3 RX1day-above threshold (",
    round(thr$threshold, 1),
    " mm, ",
    round(thr$proportion * 100, 2),
    "%)"
  )
}

# Time series of annual RX1day --------------------------------------------
plot_rx1day <- function(df, region_scenario, period_title, y_limits, thr_list, show_legend = TRUE) {
  thr <- thr_list[[region_scenario]]
  label <- threshold_label(thr)

  ggplot(df, aes(x = Year, y = RX1day)) +
    geom_line(colour = "black", linewidth = 0.35) +
    geom_hline(aes(yintercept = thr$threshold, colour = "Threshold"), linetype = "dashed", size = 1.1) +
    scale_y_continuous(limits = y_limits) +
    scale_colour_manual(values = c("Threshold" = box_colour), labels = c("Threshold" = label)) +
    labs(
      title = period_title,
      x = "Year",
      y = "RX1day (mm)",
      colour = "Threshold"
    ) +
    theme_thesis +
    theme(
      legend.position = if (show_legend) "bottom" else "none"
    )
}

# Count exceedances per year ----------------------------------------------
count_exceedances_per_year <- function(df, threshold) {
  daily_df <- df %>% select(1:360)
  apply(daily_df, 1, function(x) sum(x > threshold, na.rm = TRUE))
}

# Histogram (single threshold only) ---------------------------------------
build_hist_df <- function(region_name, thr_list, df_CD, df_FP) {
  thr_CD <- thr_list[[paste0(region_name, "_CD")]]$threshold
  thr_FP <- thr_list[[paste0(region_name, "_FP")]]$threshold

  hist_df <- bind_rows(
    data.frame(days = count_exceedances_per_year(df_CD, thr_CD), Period = "Current Day"),
    data.frame(days = count_exceedances_per_year(df_FP, thr_FP), Period = "Future Projection")
  )

  hist_df %>%
    count(Period, days) %>%
    group_by(Period) %>%
    mutate(prop_years = n / sum(n)) %>%
    ungroup() %>%
    mutate(Region = tools::toTitleCase(region_name))
}

plot_hist_exceedances <- function(hist_df_prop) {
  region_title <- unique(hist_df_prop$Region)

  ggplot(hist_df_prop, aes(x = days, y = prop_years)) +
    geom_col(width = 0.9, fill = box_colour) +
    facet_grid(. ~ Period) +
    labs(
      title = region_title,
      x = "Number of exceedance days per year",
      y = "Proportion of years"
    ) +
    theme_thesis
}

# Example years ------------------------------------------------------------
extract_daily_timeseries <- function(df, year) {
  df %>%
    filter(Year == year) %>%
    select(-Year, -RX1day) %>%
    pivot_longer(cols = everything(), names_to = "day", values_to = "rainfall_mm") %>%
    mutate(day = as.integer(day))
}

select_example_years <- function(df, threshold) {
  exceed <- count_exceedances_per_year(df, threshold)

  summary_df <- df %>%
    select(Year, RX1day) %>%
    mutate(exceedances = exceed)

  muted_year <- summary_df %>% filter(exceedances == 0) %>% arrange(RX1day) %>% slice(1) %>% pull(Year)
  single_year <- summary_df %>% filter(exceedances == 1) %>% arrange(desc(RX1day)) %>% slice(1) %>% pull(Year)
  high_year <- summary_df %>% arrange(desc(exceedances), desc(RX1day)) %>% slice(1) %>% pull(Year)

  list(muted = muted_year, single = single_year, high = high_year)
}

compute_y_max <- function(df, years, buffer = 10) {
  max_rx <- df %>% filter(Year %in% unlist(years)) %>% summarise(m = max(RX1day, na.rm = TRUE)) %>% pull(m)
  ceiling((max_rx + buffer) / 10) * 10
}

plot_daily_example <- function(df, year, threshold, title, y_max, y_lab = "") {
  ts_df <- extract_daily_timeseries(df, year)

  ggplot(ts_df, aes(x = day, y = rainfall_mm)) +
    geom_segment(aes(xend = day, y = 0, yend = rainfall_mm), colour = "black", alpha = 0.5) +
    geom_point(size = 0.7, colour = "black") +
    geom_hline(yintercept = threshold, linetype = "dashed", colour = box_colour, size = 1.1) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = title, x = "Day of Year", y = y_lab) +
    theme_thesis
}

plot_exceedance_examples <- function(df, region, period_label, threshold) {
  years <- select_example_years(df, threshold)
  y_max <- compute_y_max(df, years)

  p_muted <- plot_daily_example(df, years$muted, threshold, paste0("Muted Year (", years$muted, ")"), y_max, "Daily Rainfall (mm)")
  p_single <- plot_daily_example(df, years$single, threshold, paste0("Single Exceedance Year (", years$single, ")"), y_max)
  p_high <- plot_daily_example(df, years$high, threshold, paste0("High Exceedance Year (", years$high, ")"), y_max)

  (p_muted + p_single + p_high) +
    plot_annotation(title = paste(region, "–", period_label),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold")))
}

# Multi-panel RX1day vs exceedances ---------------------------------------
build_boxplot_panel_df <- function(regions_mod, thr_list) {
  bind_rows(lapply(regions_mod, function(reg) {
    bind_rows(lapply(c("CD", "FP"), function(period) {
      df <- get(paste0(reg, "_", period))
      threshold <- thr_list[[paste0(reg, "_", period)]]$threshold
      exceedances <- count_exceedances_per_year(df, threshold)

      data.frame(
        Region = region_labels[[reg]],
        Period = period_labels[[period]],
        RX1day = df$RX1day,
        exceedances = exceedances
      )
    }))
  }))
}

plot_rx1day_vs_exceedance_panel <- function(df_panel) {
  group_counts <- df_panel %>% count(Region, Period, exceedances, name = "n")

  df_box <- df_panel %>% left_join(group_counts, by = c("Region", "Period", "exceedances")) %>% filter(n >= 10)
  df_dot_only <- df_panel %>% left_join(group_counts, by = c("Region", "Period", "exceedances")) %>% filter(n < 10)

  df_outliers_upper <- df_box %>%
    group_by(Region, Period, exceedances) %>%
    mutate(upper_whisker = boxplot.stats(RX1day)$stats[5]) %>%
    ungroup() %>%
    filter(RX1day > upper_whisker)

  percent_per_bin <- df_panel %>%
    count(Region, Period, exceedances, name = "n") %>%
    group_by(Region, Period) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    ungroup()

  exceed_max <- max(df_panel$exceedances, na.rm = TRUE)
  rx1day_max <- max(df_panel$RX1day, na.rm = TRUE)

  ggplot() +
    geom_boxplot(
      data = df_box,
      aes(x = factor(exceedances), y = RX1day, group = exceedances),
      fill = box_colour,
      colour = box_colour,
      alpha = 0.5,
      outlier.shape = NA
    ) +
    geom_jitter(
      data = df_outliers_upper,
      aes(x = factor(exceedances), y = RX1day),
      width = 0.12,
      alpha = 0.35,
      colour = box_colour,
      size = 0.9
    ) +
    geom_jitter(
      data = df_dot_only,
      aes(x = factor(exceedances), y = RX1day),
      width = 0.12,
      alpha = 0.7,
      colour = box_colour,
      size = 1.5
    ) +
    geom_text(
      data = percent_per_bin,
      aes(
        x = factor(exceedances),
        y = rx1day_max * 1.02,
        label = ifelse(percent < 0.1, "<0.1%", paste0(round(percent, 1), "%"))
      ),
      inherit.aes = FALSE,
      size = 2.4,
      colour = "black"
    ) +
    facet_grid(Region ~ Period, switch = "y") +
    scale_x_discrete(limits = as.character(0:exceed_max)) +
    scale_y_continuous(limits = c(0, rx1day_max * 1.08)) +
    labs(x = "Number of exceedances per year", y = "RX1day (mm)", title = "RX1day vs Exceedances") +
    theme_thesis +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, size = 9, face = "bold"),
      strip.text.x = element_text(size = 10, face = "bold"),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    coord_cartesian(clip = "off")
}

# Save helper --------------------------------------------------------------
save_plot <- function(plot, filename, width = fig_width_full, height = fig_height_med) {
  ggsave(
    filename = file.path(model_data_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# Build all outputs --------------------------------------------------------
for (reg in regions_mod) {
  y_max <- max(get(paste0(reg, "_CD"))$RX1day, get(paste0(reg, "_FP"))$RX1day, na.rm = TRUE)
  y_limits <- c(0, y_max * 1.05)

  p_cd <- plot_rx1day(
    get(paste0(reg, "_CD")),
    paste0(reg, "_CD"),
    "Current Day",
    y_limits,
    thr_list,
    show_legend = TRUE
  )
  p_fp <- plot_rx1day(
    get(paste0(reg, "_FP")),
    paste0(reg, "_FP"),
    "Future Projection",
    y_limits,
    thr_list,
    show_legend = FALSE
  )

  p_ts_combined <- (p_cd | p_fp) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = region_labels[[reg]],
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
    ) &
    theme(legend.position = "bottom")

  save_plot(p_ts_combined, paste0(reg, "_rx1day_timeseries_combined.png"), width = fig_width_hoz_full, height = fig_height_med)

  p_ex_cd <- plot_exceedance_examples(
    get(paste0(reg, "_CD")),
    region_labels[[reg]],
    "Current Day",
    thr_list[[paste0(reg, "_CD")]]$threshold
  )
  p_ex_fp <- plot_exceedance_examples(
    get(paste0(reg, "_FP")),
    region_labels[[reg]],
    "Future Projection",
    thr_list[[paste0(reg, "_FP")]]$threshold
  )
  save_plot(p_ex_cd, paste0(reg, "_threshold_cd_examples.png"), height = fig_height_short)
  save_plot(p_ex_fp, paste0(reg, "_threshold_fp_examples.png"), height = fig_height_short)

  hist_df <- build_hist_df(reg, thr_list, get(paste0(reg, "_CD")), get(paste0(reg, "_FP")))
  save_plot(plot_hist_exceedances(hist_df), paste0(reg, "_exceedance_threshold_cd_fp_histogram.png"), height = fig_height_med)
}

panel_df <- build_boxplot_panel_df(regions_mod, thr_list)
panel_plot <- plot_rx1day_vs_exceedance_panel(panel_df)
save_plot(panel_plot, "all_regions_rx1day_exceedance_panel_boxplot.png", width = fig_width_hoz_full, height = fig_height_hoz_full)


# Summary tables -----------------------------------------------------------

annual_rain <- function(df) {
  day_cols <- as.character(1:360)
  df %>% mutate(AnnualRain = rowSums(across(all_of(day_cols)), na.rm = TRUE))
}

make_example_year_row <- function(df, exceedances, region, scenario, threshold_label, category, year) {
  n_years <- nrow(df)

  rx_stats <- df %>%
    arrange(desc(RX1day)) %>%
    mutate(
      RX1day_rank = row_number(),
      RX1day_percentile = (1 - (RX1day_rank - 1) / (n_years - 1)) * 100
    )

  yr_row <- df %>% filter(Year == year)
  yr_rx <- rx_stats %>% filter(Year == year)

  data.frame(
    Region = region,
    Scenario = scenario,
    Threshold = threshold_label,
    `Example Type` = category,
    Year = year,
    `Exceedance Days` = exceedances[as.character(year)],
    `RX1day (mm)` = yr_row$RX1day,
    `RX1day Rank` = yr_rx$RX1day_rank,
    `RX1day Percentile` = yr_rx$RX1day_percentile,
    `Annual Rainfall (mm)` = yr_row$AnnualRain,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

build_example_year_table <- function(df, exceedances, region, scenario, threshold_label, threshold_value) {
  yrs <- select_example_years(df, threshold_value)

  bind_rows(
    make_example_year_row(df, exceedances, region, scenario, threshold_label, "Muted", yrs$muted),
    make_example_year_row(df, exceedances, region, scenario, threshold_label, "Single Exceedance", yrs$single),
    make_example_year_row(df, exceedances, region, scenario, threshold_label, "High Exceedance", yrs$high)
  )
}

model_data_sets <- list(
  waikato_CD = annual_rain(waikato_CD),
  waikato_FP = annual_rain(waikato_FP),
  napier_CD = annual_rain(napier_CD),
  napier_FP = annual_rain(napier_FP),
  northland_CD = annual_rain(northland_CD),
  northland_FP = annual_rain(northland_FP),
  milford_CD = annual_rain(milford_CD),
  milford_FP = annual_rain(milford_FP)
)

# Exceedance vectors using each scenario's single threshold
exceedance_vectors <- lapply(names(model_data_sets), function(key) {
  count_exceedances_per_year(model_data_sets[[key]], thr_list[[key]]$threshold)
})
names(exceedance_vectors) <- names(model_data_sets)

# Example-year summary table and CSV
example_years_table <- bind_rows(lapply(names(model_data_sets), function(key) {
  parts <- strsplit(key, "_")[[1]]
  reg <- parts[1]
  scen <- parts[2]
  build_example_year_table(
    df = model_data_sets[[key]],
    exceedances = exceedance_vectors[[key]],
    region = region_labels[[reg]],
    scenario = period_labels[[scen]],
    threshold_label = "2/3 RX1day-above",
    threshold_value = thr_list[[key]]$threshold
  )
}))

write.csv(
  example_years_table,
  file = file.path(model_data_dir, "example_years_table_single_threshold.csv"),
  row.names = FALSE
)

# Mean exceedance comparison table and CSV
mean_exceedance_table <- bind_rows(lapply(regions_mod, function(reg) {
  data.frame(
    Region = region_labels[[reg]],
    `Mean Exceedance CD` = mean(exceedance_vectors[[paste0(reg, "_CD")]], na.rm = TRUE),
    `Mean Exceedance FP` = mean(exceedance_vectors[[paste0(reg, "_FP")]], na.rm = TRUE)
  )
})) %>%
  mutate(
    `Mean Change (FP - CD)` = `Mean Exceedance FP` - `Mean Exceedance CD`,
    `FP / CD` = ifelse(`Mean Exceedance CD` == 0, NA, `Mean Exceedance FP` / `Mean Exceedance CD`)
  )

write.csv(
  mean_exceedance_table,
  file = file.path(model_data_dir, "mean_exceedance_change_table_single_threshold.csv"),
  row.names = FALSE
)

# Cumulative proportions and probability ratio (FP/CD), using CD threshold for both periods
build_cumulative_exceed_table <- function(region_name, thr_list, df_CD, df_FP) {
  thr_value <- thr_list[[paste0(region_name, "_CD")]]$threshold

  exc_CD <- count_exceedances_per_year(df_CD, thr_value)
  exc_FP <- count_exceedances_per_year(df_FP, thr_value)

  max_k <- max(c(exc_CD, exc_FP), na.rm = TRUE)
  k_seq <- 0:max_k

  cum_CD <- data.frame(k = k_seq, prop_CD = sapply(k_seq, function(k) mean(exc_CD >= k)))
  cum_FP <- data.frame(k = k_seq, prop_FP = sapply(k_seq, function(k) mean(exc_FP >= k)))

  left_join(cum_CD, cum_FP, by = "k") %>%
    mutate(
      prop_CD = replace_na(prop_CD, 0),
      prop_FP = replace_na(prop_FP, 0),
      FP_to_CD = ifelse(prop_CD == 0, NA, prop_FP / prop_CD),
      Region = region_labels[[region_name]],
      Threshold = "2/3 RX1day-above (CD threshold)"
    )
}

cumulative_proportion_table <- bind_rows(lapply(regions_mod, function(reg) {
  build_cumulative_exceed_table(
    region_name = reg,
    thr_list = thr_list,
    df_CD = model_data_sets[[paste0(reg, "_CD")]],
    df_FP = model_data_sets[[paste0(reg, "_FP")]]
  )
})) %>%
  mutate(
    prop_CD = round(prop_CD, 4),
    prop_FP = round(prop_FP, 4),
    FP_to_CD = round(FP_to_CD, 4)
  ) %>%
  rename(
    `Exceedances ≥ k` = k,
    `Proportion CD` = prop_CD,
    `Proportion FP` = prop_FP,
    `FP / CD` = FP_to_CD
  ) %>%
  select(Region, Threshold, `Exceedances ≥ k`, `Proportion CD`, `Proportion FP`, `FP / CD`)

write.csv(
  cumulative_proportion_table,
  file = file.path(model_data_dir, "cumulative_proportion_table_single_threshold.csv"),
  row.names = FALSE
)
