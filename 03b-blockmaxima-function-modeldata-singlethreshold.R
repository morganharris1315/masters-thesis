# -------------------------------------------------------------------------
# 03b-blockmaxima-function-modeldata-singlethreshold.R
# -------------------------------------------------------------------------
# Feb 2026
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
box_colour <- "#93acff"
box_colour_dark <- "#6f8dff"
period_colours <- c("Current Day" = "#4C78A8", "Future Projection" = "#F58518")

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

theme_model_axes <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_line(colour = "black", linewidth = 0.3),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3)
)

# Time series of annual RX1day --------------------------------------------
plot_rx1day_combined <- function(region_name, y_limits, thr_list, df_CD, df_FP) {
  thr_cd <- thr_list[[paste0(region_name, "_CD")]]
  thr_fp <- thr_list[[paste0(region_name, "_FP")]]
  
  plot_df <- bind_rows(
    df_CD %>% select(Year, RX1day) %>% mutate(Period = "Current Day"),
    df_FP %>% select(Year, RX1day) %>% mutate(Period = "Future Projection")
  ) %>%
    mutate(Period = factor(Period, levels = c("Current Day", "Future Projection")))
  
  threshold_df <- data.frame(
    Period = factor(c("Current Day", "Future Projection"), levels = c("Current Day", "Future Projection")),
    threshold = c(thr_cd$threshold, thr_fp$threshold)
  )
  
  ggplot(plot_df, aes(x = Year, y = RX1day)) +
    geom_line(colour = "black", linewidth = 0.35) +
    geom_hline(
      data = threshold_df,
      aes(yintercept = threshold),
      colour = box_colour,
      linetype = "dashed",
      linewidth = 1.1,
      inherit.aes = FALSE
    ) +
    facet_grid(cols = vars(Period), scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = region_labels[[region_name]],
      x = "Year",
      y = "RX1day (mm)"
    ) +
    theme_thesis +
    theme_model_axes +
    theme(
      strip.placement = "outside",
      strip.text.x = element_text(face = "bold")
    )
}

# Count exceedances per year ----------------------------------------------
count_exceedances_per_year <- function(df, threshold) {
  daily_df <- df %>% select(1:360)
  apply(daily_df, 1, function(x) sum(x > threshold, na.rm = TRUE))
}

# Histogram (single threshold, fixed CD for CD/FP comparison) ---------------
build_hist_df <- function(region_name, thr_list, df_CD, df_FP) {
  thr_CD <- thr_list[[paste0(region_name, "_CD")]]$threshold
  
  hist_df <- bind_rows(
    data.frame(days = count_exceedances_per_year(df_CD, thr_CD), Period = "Current Day"),
    data.frame(days = count_exceedances_per_year(df_FP, thr_CD), Period = "Future Projection")
  )
  
  hist_df %>%
    count(Period, days) %>%
    complete(Period, days = 0:max(days, na.rm = TRUE), fill = list(n = 0)) %>%
    group_by(Period) %>%
    mutate(prop_years = n / sum(n)) %>%
    ungroup() %>%
    mutate(Region = tools::toTitleCase(region_name))
}

plot_hist_exceedances <- function(hist_df_prop, max_exceedance, fill_colour = box_colour) {
  region_title <- unique(hist_df_prop$Region)
  day_breaks <- 0:max_exceedance
  x_label_df <- hist_df_prop %>%
    distinct(Period) %>%
    mutate(
      days = mean(day_breaks),
      prop_years = 0,
      x_label = "Number of exceedance days per year"
    )
  
  ggplot(hist_df_prop, aes(x = days, y = prop_years)) +
    geom_col(width = 0.9, fill = fill_colour) +
    geom_text(
      data = x_label_df,
      aes(x = days, y = prop_years, label = x_label),
      vjust = 3.2,
      size = 3.2,
      inherit.aes = FALSE
    ) +
    facet_grid(. ~ Period) +
    scale_x_continuous(
      breaks = day_breaks,
      limits = c(-0.5, max_exceedance + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = region_title,
      x = NULL,
      y = "Proportion of years"
    ) +
    coord_cartesian(clip = "off") +
    theme_thesis +
    theme_model_axes +
    theme(
      panel.grid.major = element_line(colour = "grey82", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
      axis.line = element_line(colour = "black", linewidth = 0.2),
      plot.background = element_rect(colour = NA, fill = NA),
      plot.margin = margin(8, 10, 24, 8)
    )
}

plot_top10_exceedance_distribution <- function(hist_df_prop, max_exceedance) {
  region_title <- unique(hist_df_prop$Region)

  plot_df <- hist_df_prop %>%
    mutate(Period = factor(Period, levels = c("Current Day", "Future Projection")))

  n_label_df <- plot_df %>%
    distinct(Period, n_top10_years) %>%
    mutate(label = paste0("Top 10% years (n = ", n_top10_years, ")"))

  ggplot(plot_df, aes(x = days, y = prop_years, colour = Period)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(. ~ Period) +
    geom_text(
      data = n_label_df,
      aes(x = 0, y = Inf, label = label),
      hjust = 0,
      vjust = 1.2,
      size = 2.9,
      colour = "black",
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = period_colours, guide = "none") +
    scale_x_continuous(
      breaks = 0:max_exceedance,
      labels = c(as.character(0:(max_exceedance - 1)), paste0(max_exceedance, "+")),
      limits = c(0, max_exceedance),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = paste0(region_title, " – Distribution of Exceedance Days for Top 10% RX1day Years"),
      x = "Exceedance days per year (CD threshold)",
      y = "Share of top-10% years"
    ) +
    theme_thesis +
    theme_model_axes +
    theme(
      panel.grid.major = element_line(colour = "grey82", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
      axis.line = element_line(colour = "black", linewidth = 0.2),
      strip.text.x = element_text(face = "bold")
    )
}

# Histogram for top 10% annual RX1day years ---------------------------------
build_top_rx1day_hist_df <- function(region_name, period, thr_list, df, max_exceedance_bin = 10) {
  threshold <- thr_list[[paste0(region_name, "_CD")]]$threshold
  exceed_days <- count_exceedances_per_year(df, threshold)
  rx1day_90th <- as.numeric(quantile(df$RX1day, probs = 0.9, na.rm = TRUE, type = 7))

  top10_df <- df %>%
    transmute(
      Year,
      RX1day,
      exceed_days = exceed_days,
      exceedance_bin = pmin(exceed_days, max_exceedance_bin)
    ) %>%
    filter(RX1day >= rx1day_90th)

  hist_df <- top10_df %>%
    count(exceedance_bin, name = "n_years") %>%
    complete(exceedance_bin = 0:max_exceedance_bin, fill = list(n_years = 0)) %>%
    mutate(
      prop_years = n_years / sum(n_years),
      Period = period_labels[[period]],
      Region = region_labels[[region_name]],
      rx1day_90th = rx1day_90th,
      n_top10_years = nrow(top10_df)
    )

  hist_df
}

build_top_rx1day_hist_region_df <- function(region_name, thr_list, df_CD, df_FP, max_exceedance_bin = 10) {
  bind_rows(
    build_top_rx1day_hist_df(region_name, "CD", thr_list, df_CD, max_exceedance_bin),
    build_top_rx1day_hist_df(region_name, "FP", thr_list, df_FP, max_exceedance_bin)
  ) %>%
    rename(days = exceedance_bin)
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

day_cols <- as.character(1:360)

annual_rain <- function(df) {
  df %>%
    mutate(
      AnnualRain = rowSums(across(all_of(day_cols)), na.rm = TRUE)
    )
}

make_example_year_row <- function(df, exceedances, region, scenario, threshold_label, category, year) {
  if (!"AnnualRain" %in% colnames(df)) {
    df <- annual_rain(df)
  }

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

build_example_year_table <- function(df, threshold_value, region, scenario, threshold_label) {
  exceedances <- count_exceedances_per_year(df, threshold_value)
  names(exceedances) <- df$Year
  yrs <- select_example_years(df, threshold_value)

  bind_rows(
    make_example_year_row(df, exceedances, region, scenario, threshold_label, "Muted", yrs$muted),
    make_example_year_row(df, exceedances, region, scenario, threshold_label, "Single Exceedance", yrs$single),
    make_example_year_row(df, exceedances, region, scenario, threshold_label, "High Exceedance", yrs$high)
  )
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
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = "Day of Year", y = y_lab) +
    theme_thesis +
    theme_model_axes
}

plot_exceedance_examples <- function(df, region, period_label, threshold, y_max_override = NULL) {
  years <- select_example_years(df, threshold)
  y_max <- if (is.null(y_max_override)) compute_y_max(df, years) else y_max_override
  
  p_muted <- plot_daily_example(df, years$muted, threshold, paste0("Muted Year (", years$muted, ")"), y_max, "Daily Rainfall (mm)")
  p_single <- plot_daily_example(df, years$single, threshold, paste0("Single Exceedance Year (", years$single, ")"), y_max)
  p_high <- plot_daily_example(df, years$high, threshold, paste0("High Exceedance Year (", years$high, ")"), y_max)
  
  (p_muted + p_single + p_high) +
    plot_annotation(title = paste(region, "–", period_label),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold")))
}

compute_shared_example_y_max <- function(df_CD, thr_CD, df_FP, thr_FP) {
  years_cd <- select_example_years(df_CD, thr_CD)
  years_fp <- select_example_years(df_FP, thr_FP)
  
  max_cd <- compute_y_max(df_CD, years_cd)
  max_fp <- compute_y_max(df_FP, years_fp)
  
  max(max_cd, max_fp)
}

# Multi-panel RX1day vs exceedances (fixed CD threshold) --------------------
build_boxplot_panel_df <- function(regions_mod, thr_list) {
  bind_rows(lapply(regions_mod, function(reg) {
    bind_rows(lapply(c("CD", "FP"), function(period) {
      df <- get(paste0(reg, "_", period))
      threshold <- thr_list[[paste0(reg, "_CD")]]$threshold
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
      pct = 100 * n / total_years,
      pct_label = case_when(
        pct == 0 ~ "0%",
        pct < 0.1 ~ "<0.1%",
        TRUE ~ paste0(round(pct, 1), "%")
      ),
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
    facet_grid(Region ~ Period, switch = "y") +
    scale_x_discrete(limits = as.character(0:exceed_max)) +
    scale_y_continuous(limits = c(0, rx1day_max * 1.1)) +
    labs(x = "Number of exceedances per year", y = "RX1day (mm)", title = "RX1day vs Exceedances") +
    theme_thesis +
    theme(
      strip.placement = "outside",
      strip.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      strip.text.y.left = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      panel.spacing = unit(1.1, "lines"),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.35),
      axis.line = element_blank()
    )
}

build_timeseries_panel_df <- function(regions_mod, thr_list) {
  panel_region_order <- c("Milford Sound", "Napier", "Northland", "Waikato")

  bind_rows(lapply(regions_mod, function(reg) {
    bind_rows(lapply(c("CD", "FP"), function(period) {
      df <- get(paste0(reg, "_", period))
      threshold <- thr_list[[paste0(reg, "_", period)]]$threshold

      data.frame(
        Region = region_labels[[reg]],
        Period = period_labels[[period]],
        Year = df$Year,
        RX1day = df$RX1day,
        threshold = threshold
      )
    }))
  })) %>%
    mutate(
      Region = factor(Region, levels = panel_region_order),
      Period = factor(Period, levels = c("Current Day", "Future Projection"))
    )
}

plot_rx1day_timeseries_panel <- function(df_panel) {
  threshold_df <- df_panel %>%
    group_by(Region, Period) %>%
    summarise(
      threshold = first(threshold),
      x_pos = min(Year, na.rm = TRUE),
      y_pos = Inf,
      threshold_label = paste0("Threshold: ", round(first(threshold), 1), " mm"),
      .groups = "drop"
    )
  rx1day_max <- max(df_panel$RX1day, na.rm = TRUE)

  ggplot(df_panel, aes(x = Year, y = RX1day)) +
    geom_line(colour = "black", linewidth = 0.35) +
    geom_hline(
      data = threshold_df,
      aes(yintercept = threshold),
      colour = box_colour,
      linetype = "dashed",
      linewidth = 1.1,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = threshold_df,
      aes(x = x_pos, y = y_pos, label = threshold_label),
      hjust = 0,
      vjust = 1.15,
      size = 2.5,
      fontface = "bold",
      colour = box_colour,
      inherit.aes = FALSE
    ) +
    facet_grid(Region ~ Period, switch = "y", scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, rx1day_max * 1.05), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "RX1day Annual Time Series",
      x = "Year",
      y = "RX1day (mm)"
    ) +
    theme_thesis +
    theme(
      strip.placement = "outside",
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      strip.text.y.left = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      panel.spacing = unit(1.1, "lines"),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.35),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.3)
    )
}

# Save helper --------------------------------------------------------------
plot_output_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/single_threshold"
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width = fig_width_full, height = fig_height_med) {
  ggsave(
    filename = file.path(plot_output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

# Build all outputs --------------------------------------------------------
global_hist_max_exceedance <- max(unlist(lapply(regions_mod, function(reg) {
  thr_cd <- thr_list[[paste0(reg, "_CD")]]$threshold
  exc_cd <- count_exceedances_per_year(get(paste0(reg, "_CD")), thr_cd)
  exc_fp <- count_exceedances_per_year(get(paste0(reg, "_FP")), thr_cd)
  max(c(exc_cd, exc_fp), na.rm = TRUE)
})), na.rm = TRUE)

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
  
  shared_example_y_max <- compute_shared_example_y_max(
    get(paste0(reg, "_CD")),
    thr_list[[paste0(reg, "_CD")]]$threshold,
    get(paste0(reg, "_FP")),
    thr_list[[paste0(reg, "_FP")]]$threshold
  )
  
  p_ex_cd <- plot_exceedance_examples(
    get(paste0(reg, "_CD")),
    region_labels[[reg]],
    "Current Day",
    thr_list[[paste0(reg, "_CD")]]$threshold,
    y_max_override = shared_example_y_max
  )
  
  p_ex_fp <- plot_exceedance_examples(
    get(paste0(reg, "_FP")),
    region_labels[[reg]],
    "Future Projection",
    thr_list[[paste0(reg, "_FP")]]$threshold,
    y_max_override = shared_example_y_max
  )
  save_plot(p_ex_cd, paste0(reg, "_threshold_cd_examples.png"), width = fig_width_standard, height = fig_height_standard)
  save_plot(p_ex_fp, paste0(reg, "_threshold_fp_examples.png"), width = fig_width_standard, height = fig_height_standard)
  
  hist_df <- build_hist_df(reg, thr_list, get(paste0(reg, "_CD")), get(paste0(reg, "_FP")))
  save_plot(plot_hist_exceedances(hist_df, global_hist_max_exceedance), paste0(reg, "_exceedance_threshold_cd_fp_histogram.png"), width = fig_width_standard, height = fig_height_standard)

  top10_hist_df <- build_top_rx1day_hist_region_df(
    reg,
    thr_list,
    get(paste0(reg, "_CD")),
    get(paste0(reg, "_FP")),
    max_exceedance_bin = 10
  )
  save_plot(
    plot_top10_exceedance_distribution(top10_hist_df, max_exceedance = 10),
    paste0(reg, "_top10_rx1day_exceedance_bins_histogram.png"),
    width = fig_width_standard,
    height = fig_height_standard
  )
}

panel_df <- build_boxplot_panel_df(regions_mod, thr_list)
panel_plot <- plot_rx1day_vs_exceedance_panel(panel_df)
save_plot(panel_plot, "all_regions_rx1day_exceedance_panel_boxplot.png", width = fig_width_full, height = fig_height_tall)

timeseries_panel_df <- build_timeseries_panel_df(regions_mod, thr_list)
timeseries_panel_plot <- plot_rx1day_timeseries_panel(timeseries_panel_df)
save_plot(timeseries_panel_plot, "all_regions_rx1day_timeseries_panel.png", width = fig_width_full, height = fig_height_tall)


# Threshold summary table -------------------------------------------------

threshold_summary_table <- bind_rows(lapply(names(thr_list), function(name) {
  
  parts <- strsplit(name, "_")[[1]]
  
  data.frame(
    Region = region_labels[[parts[1]]],
    Scenario = period_labels[[parts[2]]],
    `Threshold (mm)` = round(thr_list[[name]]$threshold, 1),
    `Percentage of Daily Rainfall above Threshold` = round(thr_list[[name]]$proportion * 100, 2)
  )
}))

write.csv(
  threshold_summary_table,
  file.path(plot_output_dir, "rx1day_threshold_levels_daily_rainfall_above_threshold_percent.csv"),
  row.names = FALSE
)


# Summary table: change in mean exceedances (fixed CD threshold)

build_fixedCD_summary_table <- function(hist_df_prop) {
  hist_df_prop %>%
    group_by(Region, Period) %>%
    summarise(
      Mean_exceedances = weighted.mean(days, prop_years),
      .groups = "drop"
    )
}

fixedCD_changeinmean_table <- bind_rows(lapply(regions_mod, function(reg) {
  
  hist_df <- build_hist_df(
    reg,
    thr_list,
    get(paste0(reg, "_CD")),
    get(paste0(reg, "_FP"))
  )
  
  build_fixedCD_summary_table(hist_df)
})) %>%
  pivot_wider(
    names_from  = Period,
    values_from = Mean_exceedances
  ) %>%
  rename(
    `Mean Exceedance Days in Current Day` = `Current Day`,
    `Mean Exceedance Days in Future Projection` = `Future Projection`
  ) %>%
  mutate(
    `Change in Mean Exceedance Days` = `Mean Exceedance Days in Future Projection` - `Mean Exceedance Days in Current Day`
  )

fixedCD_changeinmean_table

write.csv(
  fixedCD_changeinmean_table,
  file.path(plot_output_dir, "mean_exceedance_change_fixedCD.csv"),
  row.names = FALSE
)

calc_cumulative_props <- function(exceed_vec) {
  max_k <- max(exceed_vec, na.rm = TRUE)
  n_years <- length(exceed_vec)
  
  data.frame(
    k = 0:max_k,
    prop_ge_k = sapply(
      0:max_k,
      function(k) sum(exceed_vec >= k, na.rm = TRUE) / n_years
    )
  )
}

build_cumulative_exceed_table_fixedCD <- function(region_name, thr_list, df_CD, df_FP) {
  
  # Fixed CD threshold 
  thr_value <- thr_list[[paste0(region_name, "_CD")]]$threshold
  
  exc_CD <- count_exceedances_per_year(df_CD, thr_value)
  exc_FP <- count_exceedances_per_year(df_FP, thr_value)
  
  max_k <- max(c(exc_CD, exc_FP), na.rm = TRUE)
  k_seq <- 0:max_k
  
  cum_CD <- data.frame(
    k = k_seq,
    prop_CD = sapply(k_seq, function(k) mean(exc_CD >= k))
  )
  
  cum_FP <- data.frame(
    k = k_seq,
    prop_FP = sapply(k_seq, function(k) mean(exc_FP >= k))
  )
  
  left_join(cum_CD, cum_FP, by = "k") %>%
    mutate(
      prop_CD = replace_na(prop_CD, 0),
      prop_FP = replace_na(prop_FP, 0),
      `FP / CD` = ifelse(prop_CD == 0, NA, prop_FP / prop_CD),
      Region = region_labels[[region_name]]
    )
}

cumulative_proportion_table <- bind_rows(lapply(regions_mod, function(reg) {
  build_cumulative_exceed_table_fixedCD(
    region_name = reg,
    thr_list = thr_list,
    df_CD = get(paste0(reg, "_CD")),
    df_FP = get(paste0(reg, "_FP"))
  )
})) %>%
  mutate(
    `Proportion in Current Day` = round(prop_CD, 4),
    `Proportion in Future Projection` = round(prop_FP, 4),
    `Probability Ratio (Future Projection/ Current Day)` = round(`FP / CD`, 4)
  ) %>%
  rename(`Exceedances ≥ k` = k) %>%
  select(
    Region,
    `Exceedances ≥ k`,
    `Proportion in Current Day`,
    `Proportion in Future Projection`,
    `Probability Ratio (Future Projection/ Current Day)`
  )

cumulative_proportion_table

write.csv(
  cumulative_proportion_table,
  file.path(plot_output_dir, "cumulative_proportion_table_fixedCD.csv"),
  row.names = FALSE
)


# Summary table: exceedance bins for top 10% RX1day years -------------------
top10_rx1day_exceedance_summary <- bind_rows(lapply(regions_mod, function(reg) {
  build_top_rx1day_hist_region_df(
    region_name = reg,
    thr_list = thr_list,
    df_CD = get(paste0(reg, "_CD")),
    df_FP = get(paste0(reg, "_FP")),
    max_exceedance_bin = 10
  )
})) %>%
  transmute(
    Region,
    Scenario = Period,
    `RX1day 90th Percentile (mm)` = round(rx1day_90th, 2),
    `Top 10% Years (count)` = n_top10_years,
    `Exceedance Days Bin` = days,
    `Number of Years in Bin` = n_years,
    `Proportion of Top 10% Years` = round(prop_years, 4)
  )

top10_rx1day_exceedance_summary

write.csv(
  top10_rx1day_exceedance_summary,
  file.path(plot_output_dir, "top10_rx1day_exceedance_bins_summary.csv"),
  row.names = FALSE
)


# Example years table ------------------------------------------------------
table_example_years_singlethreshold <- bind_rows(lapply(regions_mod, function(reg) {
  bind_rows(lapply(c("CD", "FP"), function(period) {
    df <- get(paste0(reg, "_", period)) %>% annual_rain()
    threshold_value <- thr_list[[paste0(reg, "_", period)]]$threshold

    build_example_year_table(
      df = df,
      threshold_value = threshold_value,
      region = region_labels[[reg]],
      scenario = period_labels[[period]],
      threshold_label = "Single threshold"
    )
  }))
})) %>%
  arrange(Region, Scenario, `Example Type`)

table_example_years_singlethreshold

write.csv(
  table_example_years_singlethreshold,
  file.path(plot_output_dir, "example_years_table_singlethreshold.csv"),
  row.names = FALSE
)


# Distribution plot  ------------------------------------------------------

build_rx1day_density_df <- function(region_name, period, thr_list, k = 4) {
  
  df <- get(paste0(region_name, "_", period))
  threshold <- thr_list[[paste0(region_name, "_CD")]]$threshold
  
  exceedances <- count_exceedances_per_year(df, threshold)
  
  # Full dataset: All years
  all_years_df <- data.frame(
    RX1day = df$RX1day,
    Group = "All years"
  )
  
  # Subset for years with >= k exceedances
  k_exceed_df <- data.frame(
    RX1day = df$RX1day[exceedances >= k],
    Group = paste0("≥ ", k, " exceedance days")
  )
  
  # Combine
  combined_df <- bind_rows(all_years_df, k_exceed_df) %>%
    mutate(
      Group = factor(Group, levels = c("All years", paste0("≥ ", k, " exceedance days"))),
      Region = region_labels[[region_name]],
      Period = period_labels[[period]]
    )
  
  return(combined_df)
}

plot_rx1day_density <- function(df_density, k = 4) {
  
  group_levels <- levels(df_density$Group)
  
  colour_vals <- setNames(
    c("darkgrey", "#93acff"),
    group_levels
  )
  
  # Count total number of RX1day points in all years
  n_all_years <- sum(df_density$Group == "All years")
  
  ggplot() +
    # All years
    geom_density(
      data = df_density %>% filter(Group == "All years"),
      aes(x = RX1day, colour = Group, fill = Group),
      alpha = 0.35, linewidth = 0.9, adjust = 1.1
    ) +
    # ≥ k exceedance days, scaled to proportion of all years
    geom_density(
      data = df_density %>% filter(Group != "All years"),
      aes(
        x = RX1day,
        colour = Group,
        fill = Group,
        y = ..density.. * (sum(df_density$Group != "All years") / n_all_years)
      ),
      alpha = 0.35, linewidth = 0.9, adjust = 1.1
    ) +
    scale_colour_manual(values = colour_vals) +
    scale_fill_manual(values = colour_vals) +
    labs(
      x = "RX1day (mm)",
      y = "Proportion",
      title = unique(df_density$Region),
      subtitle = unique(df_density$Period)
    ) +
    coord_cartesian(xlim = c(0, max_rx1day_all_regions)) +
    theme_thesis +
    theme_model_axes +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = c(0.78, 0.82),
      legend.justification = c(1, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.height = unit(1.1, "lines"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

k_exceed <- 4

max_rx1day_all_regions <- max(
  unlist(
    lapply(regions_mod, function(reg) {
      c(
        get(paste0(reg, "_CD"))$RX1day,
        get(paste0(reg, "_FP"))$RX1day
      )
    })
  ),
  na.rm = TRUE
)

for (reg in regions_mod) {
  for (per in c("CD", "FP")) {
    
    df_density <- build_rx1day_density_df(
      region_name = reg,
      period = per,
      thr_list = thr_list,
      k = k_exceed
    )
    
    p_density <- plot_rx1day_density(df_density, k = k_exceed)
    
    save_plot(
      p_density,
      paste0(reg, "_", tolower(per), "_rx1day_density_ge", k_exceed, "_exceedances.png"),
      height = fig_height_short
    )
  }
}

build_rx1day_density_df_all_regions <- function(period, regions_mod, thr_list, k = 4) {
  bind_rows(lapply(regions_mod, function(reg) {
    df <- get(paste0(reg, "_", period))
    threshold <- thr_list[[paste0(reg, "_CD")]]$threshold
    exceedances <- count_exceedances_per_year(df, threshold)
    
    bind_rows(
      data.frame(RX1day = df$RX1day, Group = "All years"),
      data.frame(
        RX1day = df$RX1day[exceedances >= k],
        Group = paste0("≥ ", k, " exceedance days")
      )
    )
  })) %>%
    mutate(
      Group = factor(Group, levels = c("All years", paste0("≥ ", k, " exceedance days"))),
      Region = "All Regions",
      Period = period_labels[[period]]
    )
}

for (per in c("CD", "FP")) {
  df_density_all_regions <- build_rx1day_density_df_all_regions(
    period = per,
    regions_mod = regions_mod,
    thr_list = thr_list,
    k = k_exceed
  )
  
  p_density_all_regions <- plot_rx1day_density(df_density_all_regions, k = k_exceed)
  
  save_plot(
    p_density_all_regions,
    paste0("all_regions_", tolower(per), "_rx1day_density_ge", k_exceed, "_exceedances.png"),
    height = fig_height_short
  )
}
