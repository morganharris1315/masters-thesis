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
box_colour <- "#6A5ACD"

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
plot_rx1day_combined <- function(region_name, y_limits, thr_list, df_CD, df_FP) {
  thr_cd <- thr_list[[paste0(region_name, "_CD")]]
  thr_fp <- thr_list[[paste0(region_name, "_FP")]]
  label <- c(
    paste0("Current Day: ", threshold_label(thr_cd)),
    paste0("Future Projection: ", threshold_label(thr_fp))
  )

  plot_df <- bind_rows(
    df_CD %>% select(Year, RX1day) %>% mutate(Period = "Current Day"),
    df_FP %>% select(Year, RX1day) %>% mutate(Period = "Future Prediction")
  ) %>%
    mutate(Period = factor(Period, levels = c("Current Day", "Future Prediction")))

  threshold_df <- data.frame(
    Period = factor(c("Current Day", "Future Prediction"), levels = c("Current Day", "Future Prediction")),
    threshold = c(thr_cd$threshold, thr_fp$threshold)
  )

  ggplot(plot_df, aes(x = Year, y = RX1day)) +
    geom_line(colour = "black", linewidth = 0.35) +
    geom_hline(
      data = threshold_df,
      aes(yintercept = threshold, colour = "Threshold"),
      linetype = "dashed",
      linewidth = 1.1,
      inherit.aes = FALSE
    ) +
    facet_grid(rows = vars(Period)) +
    scale_y_continuous(limits = y_limits) +
    scale_colour_manual(values = c("Threshold" = box_colour), labels = c("Threshold" = paste(label, collapse = "\n"))) +
    labs(
      title = region_labels[[region_name]],
      x = "Year",
      y = "RX1day (mm)",
      colour = "Threshold"
    ) +
    theme_thesis +
    theme(
      strip.text.y = element_text(angle = 0)
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
    data.frame(days = count_exceedances_per_year(df_FP, thr_FP), Period = "Future Prediction")
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
  period_totals <- df_panel %>% count(Region, Period, name = "total_years")
  
  df_box <- df_panel %>% left_join(group_counts, by = c("Region", "Period", "exceedances")) %>% filter(n >= 10)
  df_dot_only <- df_panel %>% left_join(group_counts, by = c("Region", "Period", "exceedances")) %>% filter(n < 10)

  box_stats <- df_box %>%
    group_by(Region, Period, exceedances) %>%
    summarise(
      q1 = quantile(RX1day, 0.25, na.rm = TRUE),
      q3 = quantile(RX1day, 0.75, na.rm = TRUE),
      iqr = IQR(RX1day, na.rm = TRUE),
      upper_whisker = max(RX1day[RX1day <= q3 + 1.5 * iqr], na.rm = TRUE),
      .groups = "drop"
    )

  df_box_dots_above <- df_box %>%
    left_join(box_stats, by = c("Region", "Period", "exceedances")) %>%
    filter(RX1day > upper_whisker)
  
  exceed_max <- max(df_panel$exceedances, na.rm = TRUE)
  rx1day_max <- max(df_panel$RX1day, na.rm = TRUE)

  pct_labels <- group_counts %>%
    left_join(period_totals, by = c("Region", "Period")) %>%
    mutate(
      pct_label = paste0(round(100 * n / total_years, 1), "%"),
      y = rx1day_max * 1.03
    )
  
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
      data = df_box_dots_above,
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
      data = pct_labels,
      aes(x = factor(exceedances), y = y, label = pct_label),
      size = 2.5,
      vjust = 0,
      colour = "black"
    ) +
    geom_vline(
      xintercept = c(0.5, exceed_max + 0.5),
      colour = "grey45",
      linewidth = 0.35
    ) +
    geom_hline(
      yintercept = c(0, rx1day_max * 1.1),
      colour = "grey45",
      linewidth = 0.35
    ) +
    facet_grid(Region ~ Period, switch = "y") +
    scale_x_discrete(limits = as.character(0:exceed_max)) +
    scale_y_continuous(limits = c(0, rx1day_max * 1.1)) +
    labs(x = "Number of exceedances per year", y = "RX1day (mm)", title = "RX1day vs Exceedances") +
    theme_thesis +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 0),
      panel.spacing = unit(1.1, "lines"),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.35),
      axis.line = element_line(colour = "grey35", linewidth = 0.3)
    )
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
  
  p_rx1day <- plot_rx1day_combined(
    reg,
    y_limits,
    thr_list,
    get(paste0(reg, "_CD")),
    get(paste0(reg, "_FP"))
  )
  save_plot(p_rx1day, paste0(reg, "_rx1day_timeseries_combined.png"), height = fig_height_med)
  
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


# Threshold summary table -------------------------------------------------

threshold_summary_table <- bind_rows(lapply(names(thr_list), function(name) {
  
  parts <- strsplit(name, "_")[[1]]
  
  data.frame(
    Region = region_labels[[parts[1]]],
    Period = period_labels[[parts[2]]],
    Threshold_definition = "2/3 RX1day-above",
    Threshold_mm = thr_list[[name]]$threshold,
    Proportion_of_days_exceeding = thr_list[[name]]$proportion
  )
}))

write.csv(
  threshold_summary_table,
  file.path(model_data_dir, "rx1day_single_threshold_summary.csv"),
  row.names = FALSE
)


# Mean exceedances per year (CD vs FP) ------------------------------------

build_mean_exceedance_table <- function(region_name, thr_list, df_CD, df_FP) {
  
  thr_CD <- thr_list[[paste0(region_name, "_CD")]]$threshold
  thr_FP <- thr_list[[paste0(region_name, "_FP")]]$threshold
  
  mean_CD <- mean(count_exceedances_per_year(df_CD, thr_CD), na.rm = TRUE)
  mean_FP <- mean(count_exceedances_per_year(df_FP, thr_FP), na.rm = TRUE)
  
  data.frame(
    Region = region_labels[[region_name]],
    Mean_exceedances_CD = mean_CD,
    Mean_exceedances_FP = mean_FP,
    Absolute_change = mean_FP - mean_CD,
    Relative_change = mean_FP / mean_CD
  )
}

mean_exceedance_summary <- bind_rows(lapply(regions_mod, function(reg) {
  build_mean_exceedance_table(
    reg,
    thr_list,
    get(paste0(reg, "_CD")),
    get(paste0(reg, "_FP"))
  )
}))

write.csv(
  mean_exceedance_summary,
  file.path(model_data_dir, "rx1day_single_threshold_mean_exceedances.csv"),
  row.names = FALSE
)


# Cumulative exceedance probability table ---------------------------------

build_cumulative_exceedance_table <- function(region_name, thr_list, df_CD, df_FP) {
  
  thr_CD <- thr_list[[paste0(region_name, "_CD")]]$threshold
  thr_FP <- thr_list[[paste0(region_name, "_FP")]]$threshold
  
  bind_rows(
    data.frame(
      exceedances = count_exceedances_per_year(df_CD, thr_CD),
      Period = "Current Day"
    ),
    data.frame(
      exceedances = count_exceedances_per_year(df_FP, thr_FP),
      Period = "Future Projection"
    )
  ) %>%
    count(Period, exceedances) %>%
    group_by(Period) %>%
    mutate(Proportion_years = n / sum(n)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = Period,
      values_from = Proportion_years
    ) %>%
    mutate(
      Region = region_labels[[region_name]],
      FP_over_CD = `Future Projection` / `Current Day`
    ) %>%
    relocate(Region, exceedances)
}

cumulative_exceedance_table <- bind_rows(lapply(regions_mod, function(reg) {
  build_cumulative_exceedance_table(
    reg,
    thr_list,
    get(paste0(reg, "_CD")),
    get(paste0(reg, "_FP"))
  )
}))

write.csv(
  cumulative_exceedance_table,
  file.path(model_data_dir, "rx1day_single_threshold_cumulative_exceedance.csv"),
  row.names = FALSE
)
