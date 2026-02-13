
# -------------------------------------------------------------------------
# 03b-blockmaxima-function-modeldata-multithreshold.R
# -------------------------------------------------------------------------
# Jan 2026
# Creating functions to make the processing of model data easier
# -------------------------------------------------------------------------

# Read raw data ------------------------------------------------------------
waikato_CD_raw<- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/WaikatoCurrentClimate.csv")
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
    mutate(RX1day = apply(across(1:360), 1, max, na.rm = TRUE), Year = row_number())}

waikato_CD   <- prepare_model_data(waikato_CD_raw)
waikato_FP   <- prepare_model_data(waikato_FP_raw)
napier_CD    <- prepare_model_data(napier_CD_raw)
napier_FP    <- prepare_model_data(napier_FP_raw)
northland_CD <- prepare_model_data(northland_CD_raw)
northland_FP <- prepare_model_data(northland_FP_raw)
milford_CD <- prepare_model_data(milford_CD_raw)
milford_FP <- prepare_model_data(milford_FP_raw)

# Calculate thresholds ----------------------------------------------------
calculate_rx1day_thresholds <- function(df) {
  
  rx <- df$RX1day
  min_rx  <- min(rx, na.rm = TRUE)
  mean_rx <- mean(rx, na.rm = TRUE)
  max_rx  <- max(rx, na.rm = TRUE)
  mid_rx  <- (min_rx + mean_rx) / 2
  
  daily_vec <- df %>%
    select(1:360) %>%
    as.matrix() %>%
    as.numeric()
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
      min       = min_rx,
      mid       = mid_rx,
      mean      = mean_rx,
      max       = max_rx,
      onein360  = as.numeric(threshold_onein360)
    ),
    proportions = c(prop_exceed, onein360 = prop_onein360)
  )
}

thr_waikato_CD <- calculate_rx1day_thresholds(waikato_CD)
thr_waikato_FP <- calculate_rx1day_thresholds(waikato_FP)
thr_napier_CD  <- calculate_rx1day_thresholds(napier_CD)
thr_napier_FP  <- calculate_rx1day_thresholds(napier_FP)
thr_northland_CD <- calculate_rx1day_thresholds(northland_CD)
thr_northland_FP <- calculate_rx1day_thresholds(northland_FP)
thr_milford_CD<- calculate_rx1day_thresholds(milford_CD)
thr_milford_FP <- calculate_rx1day_thresholds(milford_FP)


#checking it runs correctly
thr_waikato_CD$thresholds
thr_waikato_FP$proportions

#threshold into a list
thr_list <- list(
  waikato_CD   = calculate_rx1day_thresholds(waikato_CD),
  waikato_FP   = calculate_rx1day_thresholds(waikato_FP),
  napier_CD    = calculate_rx1day_thresholds(napier_CD),
  napier_FP    = calculate_rx1day_thresholds(napier_FP),
  northland_CD = calculate_rx1day_thresholds(northland_CD),
  northland_FP = calculate_rx1day_thresholds(northland_FP),
  milford_CD = calculate_rx1day_thresholds(milford_CD),
  milford_FP = calculate_rx1day_thresholds(milford_FP)
)

make_labels <- function(region_scenario, thr_list) {
  
  thr <- thr_list[[region_scenario]]  # region_scenario example: "waikato_CD"
  
  labels <- c(
    "Min RX1day" = paste0(
      "Min RX1day (", round(thr$thresholds["min"], 1), " mm, ",
      round(thr$proportions["min"] * 100, 2), "%)"
    ),
    "Mid RX1day" = paste0(
      "Mid RX1day (", round(thr$thresholds["mid"], 1), " mm, ",
      round(thr$proportions["mid"] * 100, 2), "%)"
    ),
    "Mean RX1day" = paste0(
      "Mean RX1day (", round(thr$thresholds["mean"], 1), " mm, ",
      round(thr$proportions["mean"] * 100, 2), "%)"
    ),
    "1 in 360 day threshold" = paste0(
      "1 in 360 day (", round(thr$thresholds["onein360"], 1), " mm, ",
      round(thr$proportions["onein360"] * 100, 2), "%)"
    ),
    "Max RX1day" = paste0(
      "Max RX1day (", round(thr$thresholds["max"], 1), " mm, ",
      round(thr$proportions["max"] * 100, 2), "%)"
    )
  )
  
  return(labels)
}

labels_waikato_CD   <- make_labels("waikato_CD", thr_list)
labels_waikato_FP   <- make_labels("waikato_FP", thr_list)
labels_napier_CD    <- make_labels("napier_CD", thr_list)
labels_napier_FP    <- make_labels("napier_FP", thr_list)
labels_northland_CD <- make_labels("northland_CD", thr_list)
labels_northland_FP <- make_labels("northland_FP", thr_list)
labels_milford_CD <- make_labels("milford_CD", thr_list)
labels_milford_FP <- make_labels("milford_FP", thr_list)


# RX1day timeseries plot --------------------------------------------------
plot_rx1day <- function(df, region_scenario, labels, y_limits, thr_list) {
  
  thr <- thr_list[[region_scenario]]
  
  parts <- strsplit(region_scenario, "_")[[1]]
  region <- tools::toTitleCase(parts[1])  # Capitalize first letter
  scenario <- ifelse(parts[2] == "CD", "Current Day", "Future Projections")
  plot_title <- paste(region, "-", scenario)
  
  ggplot(df, aes(x = Year, y = RX1day)) +
    geom_line() +
    theme_thesis +
    labs(
      title = plot_title,
      x = "Year",
      y = "RX1day (mm)",
      colour = "Threshold"
    ) +
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
    )
}

y_max_waikato <- max(waikato_CD$RX1day, waikato_FP$RX1day, na.rm = TRUE)
y_limits_waikato <- c(0, y_max_waikato * 1.05)  # 5% headroom

y_max_napier <- max(napier_CD$RX1day, napier_FP$RX1day, na.rm = TRUE)
y_limits_napier <- c(0, y_max_napier * 1.05)

y_max_northland <- max(northland_CD$RX1day, northland_FP$RX1day, na.rm = TRUE)
y_limits_northland <- c(0, y_max_northland * 1.05)

y_max_milford <- max(milford_CD$RX1day, milford_FP$RX1day, na.rm = TRUE)
y_limits_milford <- c(0, y_max_milford * 1.05)

p_RX1day_timseries_waikato_CD <- plot_rx1day(waikato_CD, "waikato_CD", labels_waikato_CD, y_limits_waikato, thr_list)
p_RX1day_timseries_waikato_FP <- plot_rx1day(waikato_FP, "waikato_FP", labels_waikato_FP, y_limits_waikato, thr_list)

p_RX1day_timseries_napier_CD <- plot_rx1day(napier_CD, "napier_CD", labels_napier_CD, y_limits_napier, thr_list)
p_RX1day_timseries_napier_FP <- plot_rx1day(napier_FP, "napier_FP", labels_napier_FP, y_limits_napier, thr_list)

p_RX1day_timseries_northland_CD <- plot_rx1day(northland_CD, "northland_CD", labels_northland_CD, y_limits_northland, thr_list)
p_RX1day_timseries_northland_FP <- plot_rx1day(northland_FP, "northland_FP", labels_northland_FP, y_limits_northland, thr_list)

p_RX1day_timseries_milford_CD <- plot_rx1day(milford_CD, "milford_CD", labels_milford_CD, y_limits_milford, thr_list)
p_RX1day_timseries_milford_FP <- plot_rx1day(milford_FP, "milford_FP", labels_milford_FP, y_limits_milford, thr_list)

#checking plots
p_RX1day_timseries_waikato_CD

p_RX1day_timseries_waikato_combined <- p_RX1day_timseries_waikato_CD / p_RX1day_timseries_waikato_FP
p_RX1day_timseries_napier_combined <- p_RX1day_timseries_napier_CD / p_RX1day_timseries_napier_FP
p_RX1day_timseries_northland_combined <- p_RX1day_timseries_northland_CD / p_RX1day_timseries_northland_FP
p_RX1day_timseries_milford_combined <- p_RX1day_timseries_milford_CD / p_RX1day_timseries_milford_FP


# Count exceedances per year ----------------------------------------------
count_exceedances_per_year <- function(df, threshold) {
  daily_df <- df %>% select(1:360)
  
  exc <- apply(daily_df, 1, function(x) {
    sum(x > threshold, na.rm = TRUE)
  })
  
  names(exc) <- df$Year 
  return(exc)
}


build_hist_df <- function(region_name, thr_list, df_CD, df_FP) {
  
  thr_CD <- thr_list[[paste0(region_name, "_CD")]]$thresholds
  thr_FP <- thr_list[[paste0(region_name, "_FP")]]$thresholds
  
  hist_df <- bind_rows(
    
    # Current Day
    data.frame(
      days = count_exceedances_per_year(df_CD, thr_CD["mid"]),
      Threshold = "> Mid RX1day",
      Period = "Current Day"
    ),
    data.frame(
      days = count_exceedances_per_year(df_CD, thr_CD["onein360"]),
      Threshold = "> 1-in-360 threshold",
      Period = "Current Day"
    ),
    data.frame(
      days = count_exceedances_per_year(df_CD, thr_CD["mean"]),
      Threshold = "> Mean RX1day",
      Period = "Current Day"
    ),
    
    # Future Projection
    data.frame(
      days = count_exceedances_per_year(df_FP, thr_FP["mid"]),
      Threshold = "> Mid RX1day",
      Period = "Future Projection"
    ),
    data.frame(
      days = count_exceedances_per_year(df_FP, thr_FP["onein360"]),
      Threshold = "> 1-in-360 threshold",
      Period = "Future Projection"
    ),
    data.frame(
      days = count_exceedances_per_year(df_FP, thr_FP["mean"]),
      Threshold = "> Mean RX1day",
      Period = "Future Projection"
    )
  )
  
  # Factor ordering
  hist_df$Threshold <- factor(
    hist_df$Threshold,
    levels = c(
      "> Mid RX1day",
      "> 1-in-360 threshold",
      "> Mean RX1day"
    )
  )
  
  # Convert to proportions per facet
  hist_df_prop <- hist_df %>%
    count(Period, Threshold, days) %>%
    group_by(Period, Threshold) %>%
    mutate(prop_years = n / sum(n)) %>%
    ungroup()
  
  # Store region name for plotting
  hist_df_prop$Region <- tools::toTitleCase(region_name)
  
  hist_df_prop
}

plot_hist_exceedances <- function(hist_df_prop) {
  
  region_title <- unique(hist_df_prop$Region)
  
  ggplot(hist_df_prop, aes(x = days, y = prop_years, fill = Threshold)) +
    geom_col(width = 0.9) +
    facet_grid(Period ~ Threshold, switch = "y") +
    scale_x_continuous(
      limits = c(-0.5, 12.5),
      breaks = 0:12
    ) +
    scale_y_continuous(
      limits = c(0, 0.6),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = c(
        "> Mid RX1day" = "#0072B2",
        "> 1-in-360 threshold" = "#E69F00",
        "> Mean RX1day" = "#009E73"
      ),
      guide = "none"
    ) +
    labs(
      title = region_title,
      x = "Number of exceedance days per year",
      y = "Proportion of years"
    ) +
    theme_thesis +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90),
      panel.spacing = unit(0.8, "lines"),
      plot.title = element_text(hjust = 0.5)
    )
}

hist_df_waikato   <- build_hist_df("waikato", thr_list, waikato_CD, waikato_FP)
hist_plot_waikato <- plot_hist_exceedances(hist_df_waikato)

hist_df_napier    <- build_hist_df("napier", thr_list, napier_CD, napier_FP)
hist_plot_napier  <- plot_hist_exceedances(hist_df_napier)

hist_df_northland <- build_hist_df("northland", thr_list, northland_CD, northland_FP)
hist_plot_northland <- plot_hist_exceedances(hist_df_northland)

hist_df_milford <- build_hist_df("milford", thr_list, milford_CD, milford_FP)
hist_plot_milford <- plot_hist_exceedances(hist_df_milford)

# example years -----------------------------------------------------------
extract_daily_timeseries <- function(df, year) {
  
  df %>%
    filter(Year == year) %>%
    select(-Year, -RX1day) %>%
    pivot_longer(
      cols = everything(),
      names_to = "day",
      values_to = "rainfall_mm"
    ) %>%
    mutate(day = as.integer(day))
}

select_example_years <- function(df, threshold) {
  
  exceed <- count_exceedances_per_year(df, threshold)
  
  summary_df <- df %>%
    select(Year, RX1day) %>%
    mutate(exceedances = exceed)
  
  muted_year <- summary_df %>%
    filter(exceedances == 0) %>%
    arrange(RX1day) %>%
    slice(1) %>%
    pull(Year)
  
  single_year <- summary_df %>%
    filter(exceedances == 1) %>%
    arrange(desc(RX1day)) %>%
    slice(1) %>%
    pull(Year)
  
  high_year <- summary_df %>%
    arrange(desc(exceedances), desc(RX1day)) %>%
    slice(1) %>%
    pull(Year)
  
  list(
    muted  = muted_year,
    single = single_year,
    high   = high_year
  )
}
compute_y_max <- function(df, years, buffer = 10) {
  
  max_rx <- df %>%
    filter(Year %in% unlist(years)) %>%
    summarise(max_RX1day = max(RX1day, na.rm = TRUE)) %>%
    pull(max_RX1day)
  
  y_max_raw <- max_rx + buffer
  
  ceiling(y_max_raw / 10) * 10
}
plot_daily_example <- function(df, year, threshold, title, y_max, colour, y_lab = "") {
  
  ts_df <- extract_daily_timeseries(df, year)
  
  ggplot(ts_df, aes(x = day, y = rainfall_mm)) +
    geom_line() +
    geom_hline(
      yintercept = threshold,
      linetype = "dashed",
      colour = colour,
      size = 1.1
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, by = 50)
    ) +
    labs(
      title = title,
      x = "Day of Year",
      y = y_lab
    ) +
    theme_thesis +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}
plot_exceedance_examples <- function(df, region, period_label, threshold, threshold_label, colour) {
  
  years <- select_example_years(df, threshold)
  y_max <- compute_y_max(df, years)
  
  p_muted <- plot_daily_example(
    df, years$muted, threshold,
    paste("Muted Year (", years$muted, ")", sep = ""),
    y_max, colour,
    y_lab = "Daily Rainfall (mm)"
  )
  
  p_single <- plot_daily_example(
    df, years$single, threshold,
    paste("Single Exceedance Year (", years$single, ")", sep = ""),
    y_max, colour
  )
  
  p_high <- plot_daily_example(
    df, years$high, threshold,
    paste("High Exceedance Year (", years$high, ")", sep = ""),
    y_max, colour
  )
  
  (p_muted + p_single + p_high) +
    plot_annotation(
      title = paste(region, "–", period_label, "–", threshold_label),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 9, face = "bold")
      )
    )
}

#Checking years
select_example_years(waikato_CD, thr_list[["waikato_CD"]]$thresholds["mid"])

p_exampleyears_waikato_CD_mid <- plot_exceedance_examples(waikato_CD, "Waikato", "Current Day",
                                                          thr_list[["waikato_CD"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_waikato_CD_360 <- plot_exceedance_examples(waikato_CD, "Waikato", "Current Day",
                                                          thr_list[["waikato_CD"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

p_exampleyears_waikato_FP_mid <- plot_exceedance_examples(waikato_FP, "Waikato", "Future Projection",
                                                          thr_list[["waikato_FP"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_waikato_FP_360 <- plot_exceedance_examples(waikato_FP, "Waikato", "Future Projection",
                                                          thr_list[["waikato_FP"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

p_exampleyears_napier_CD_mid <- plot_exceedance_examples(napier_CD, "Napier", "Current Day",
                                                         thr_list[["napier_CD"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_napier_CD_360 <- plot_exceedance_examples(napier_CD, "Napier", "Current Day",
                                                         thr_list[["napier_CD"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

p_exampleyears_napier_FP_mid <- plot_exceedance_examples(napier_FP, "Napier", "Future Projection",
                                                         thr_list[["napier_FP"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_napier_FP_360 <- plot_exceedance_examples(napier_FP, "Napier", "Future Projection",
                                                         thr_list[["napier_FP"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

p_exampleyears_northland_CD_mid <- plot_exceedance_examples(northland_CD, "Northland", "Current Day",
                                                            thr_list[["northland_CD"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_northland_CD_360 <- plot_exceedance_examples(northland_CD, "Northland", "Current Day",
                                                            thr_list[["northland_CD"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

p_exampleyears_northland_FP_mid <- plot_exceedance_examples(northland_FP, "Northland", "Future Projection",
                                                            thr_list[["northland_FP"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_northland_FP_360 <- plot_exceedance_examples(northland_FP, "Northland", "Future Projection",
                                                            thr_list[["northland_FP"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

p_exampleyears_milford_CD_mid <- plot_exceedance_examples(milford_CD, "Milford Sound", "Current Day",
                                                          thr_list[["milford_CD"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_milford_CD_360 <- plot_exceedance_examples(milford_CD, "Milford Sound", "Current Day",
                                                          thr_list[["milford_CD"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

p_exampleyears_milford_FP_mid <- plot_exceedance_examples(milford_FP, "Milford Sound", "Future Projection",
                                                          thr_list[["milford_FP"]]$thresholds["mid"], "Mid RX1day", "#0072B2")

p_exampleyears_milford_FP_360 <- plot_exceedance_examples(milford_FP, "Milford Sound", "Future Projection",
                                                          thr_list[["milford_FP"]]$thresholds["onein360"], "1-in-360 Threshold", "#E69F00")

# Comparing RX1day and exceedances ----------------------------------------
waikato_cd_mid <- count_exceedances_per_year(waikato_CD, thr_list[["waikato_CD"]]$thresholds["mid"])
waikato_cd_360 <- count_exceedances_per_year(waikato_CD, thr_list[["waikato_CD"]]$thresholds["onein360"])
waikato_fp_mid <- count_exceedances_per_year(waikato_FP, thr_list[["waikato_FP"]]$thresholds["mid"])
waikato_fp_360 <- count_exceedances_per_year(waikato_FP, thr_list[["waikato_FP"]]$thresholds["onein360"])

napier_cd_mid <- count_exceedances_per_year(napier_CD, thr_list[["napier_CD"]]$thresholds["mid"])
napier_cd_360 <- count_exceedances_per_year(napier_CD, thr_list[["napier_CD"]]$thresholds["onein360"])
napier_fp_mid <- count_exceedances_per_year(napier_FP, thr_list[["napier_FP"]]$thresholds["mid"])
napier_fp_360 <- count_exceedances_per_year(napier_FP, thr_list[["napier_FP"]]$thresholds["onein360"])

northland_cd_mid <- count_exceedances_per_year(northland_CD, thr_list[["northland_CD"]]$thresholds["mid"])
northland_cd_360 <- count_exceedances_per_year(northland_CD, thr_list[["northland_CD"]]$thresholds["onein360"])
northland_fp_mid <- count_exceedances_per_year(northland_FP, thr_list[["northland_FP"]]$thresholds["mid"])
northland_fp_360 <- count_exceedances_per_year(northland_FP, thr_list[["northland_FP"]]$thresholds["onein360"])

milford_cd_mid <- count_exceedances_per_year(milford_CD, thr_list[["milford_CD"]]$thresholds["mid"])
milford_cd_360 <- count_exceedances_per_year(milford_CD, thr_list[["milford_CD"]]$thresholds["onein360"])
milford_fp_mid <- count_exceedances_per_year(milford_FP, thr_list[["milford_FP"]]$thresholds["mid"])
milford_fp_360 <- count_exceedances_per_year(milford_FP, thr_list[["milford_FP"]]$thresholds["onein360"])

waikato_rx1day_max <- max(c(waikato_CD$RX1day, waikato_FP$RX1day), na.rm = TRUE) + 10
waikato_exceed_max <- max(waikato_cd_mid, waikato_cd_360, waikato_fp_mid, waikato_fp_360, na.rm = TRUE)
napier_rx1day_max <- max(c(napier_CD$RX1day, napier_FP$RX1day), na.rm = TRUE) + 10
napier_exceed_max <- max(napier_cd_mid, napier_cd_360, napier_fp_mid, napier_fp_360, na.rm = TRUE)
northland_rx1day_max <- max(c(northland_CD$RX1day, northland_FP$RX1day), na.rm = TRUE) + 10
northland_exceed_max <- max(northland_cd_mid, northland_cd_360, northland_fp_mid, northland_fp_360, na.rm = TRUE)
milford_rx1day_max <- max(c(milford_CD$RX1day, milford_FP$RX1day), na.rm = TRUE) + 10
milford_exceed_max <- max(milford_cd_mid, milford_cd_360, milford_fp_mid, milford_fp_360, na.rm = TRUE)

plot_rx1day_vs_exceedance <- function(data, exceedances, title_text) {
  region_name <- deparse(substitute(data))
  region <- sub("_.*", "", region_name)
  rx1day_max <- get(paste0(region, "_rx1day_max"))
  exceed_max <- get(paste0(region, "_exceed_max"))
  data %>%
    mutate(exceedances = exceedances) %>%
    select(Year, RX1day, exceedances) %>%
    ggplot(aes(x = factor(exceedances), y = RX1day, fill = factor(exceedances))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4) +
    theme_thesis +
    labs(x = "Number of exceedances per year", y = "RX1day (mm)", title = title_text) +
    guides(fill = "none") +
    scale_y_continuous(limits = c(0, rx1day_max), breaks = seq(0, rx1day_max, by = 100)) +
    scale_x_discrete(limits = as.character(0:exceed_max))
}

boxplot_cd_mid_waikato <- plot_rx1day_vs_exceedance(waikato_CD, waikato_cd_mid, "Mid RX1day - Current Day")
boxplot_cd_360_waikato <- plot_rx1day_vs_exceedance(waikato_CD, waikato_cd_360, "1-in-360 - Current Day")
boxplot_fp_mid_waikato <- plot_rx1day_vs_exceedance(waikato_FP, waikato_fp_mid, "Mid RX1day - Future Projection")
boxplot_fp_360_waikato <- plot_rx1day_vs_exceedance(waikato_FP, waikato_fp_360, "1-in-360 - Future Projection")

boxplot_cd_mid_napier <- plot_rx1day_vs_exceedance(napier_CD, napier_cd_mid, "Mid RX1day - Current Day")
boxplot_cd_360_napier <- plot_rx1day_vs_exceedance(napier_CD, napier_cd_360, "1-in-360 - Current Day")
boxplot_fp_mid_napier <- plot_rx1day_vs_exceedance(napier_FP, napier_fp_mid, "Mid RX1day - Future Projection")
boxplot_fp_360_napier <- plot_rx1day_vs_exceedance(napier_FP, napier_fp_360, "1-in-360 - Future Projection")

boxplot_cd_mid_northland <- plot_rx1day_vs_exceedance(northland_CD, northland_cd_mid, "Mid RX1day - Current Day")
boxplot_cd_360_northland <- plot_rx1day_vs_exceedance(northland_CD, northland_cd_360, "1-in-360 - Current Day")
boxplot_fp_mid_northland <- plot_rx1day_vs_exceedance(northland_FP, northland_fp_mid, "Mid RX1day - Future Projection")
boxplot_fp_360_northland <- plot_rx1day_vs_exceedance(northland_FP, northland_fp_360, "1-in-360 - Future Projection")

boxplot_cd_mid_milford <- plot_rx1day_vs_exceedance(milford_CD, milford_cd_mid, "Mid RX1day - Current Day")
boxplot_cd_360_milford <- plot_rx1day_vs_exceedance(milford_CD, milford_cd_360, "1-in-360 - Current Day")
boxplot_fp_mid_milford <- plot_rx1day_vs_exceedance(milford_FP, milford_fp_mid, "Mid RX1day - Future Projection")
boxplot_fp_360_milford <- plot_rx1day_vs_exceedance(milford_FP, milford_fp_360, "1-in-360 - Future Projection")

# General function to save plots
model_data_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/multi_threshold"
dir.create(model_data_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width = fig_width_full, height = fig_height_med) {
  ggsave(filename = file.path(model_data_dir, filename),
         plot = plot,
         width = width,
         height = height,
         dpi = 300)
}

# List of regions
regions_mod <- c("waikato", "napier", "northland", "milford")
region_labels <- c(
  waikato = "Waikato",
  napier = "Napier",
  northland = "Northland",
  milford = "Milford Sound"
)

for (reg in regions_mod) {
  save_plot(get(paste0("p_RX1day_timseries_", reg, "_combined")), paste0(reg, "_rx1day_timeseries_combined.png"),height = fig_height_med)
  save_plot(get(paste0("p_exampleyears_", reg, "_CD_mid")), paste0(reg, "_mid_cd_examples.png"), width = fig_width_standard, height = fig_height_standard)
  save_plot(get(paste0("p_exampleyears_", reg, "_FP_mid")), paste0(reg, "_mid_fp_examples.png"), width = fig_width_standard, height = fig_height_standard)
  save_plot(get(paste0("p_exampleyears_", reg, "_CD_360")), paste0(reg, "_1in360_cd_examples.png"), width = fig_width_standard, height = fig_height_standard)
  save_plot(get(paste0("p_exampleyears_", reg, "_FP_360")), paste0(reg, "_1in360_fp_examples.png"), width = fig_width_standard, height = fig_height_standard)
  
  box_cd_mid  <- get(paste0("boxplot_cd_mid_", reg))
  box_fp_mid  <- get(paste0("boxplot_fp_mid_", reg))
  box_cd_360  <- get(paste0("boxplot_cd_360_", reg))
  box_fp_360  <- get(paste0("boxplot_fp_360_", reg))
  
  combined_boxplot <-(box_cd_mid | box_fp_mid) / (box_cd_360 | box_fp_360) +
    plot_annotation(title = paste(region_labels[reg], "- RX1day Exceedances"),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold")))
  
  save_plot(combined_boxplot, paste0(reg, "_rx1day_all_boxplots.png"), width = fig_width_full,height = fig_height_tall)
  
  save_plot(get(paste0("hist_plot_", reg)), paste0(reg, "_exceedance_allthresholds_cd_fp_histogram.png"), width = fig_width_standard, height = fig_height_standard)
}

# Summary tables ----------------------------------------------------------

#Example years
select_example_years(waikato_CD,  thr_list[["waikato_CD"]]$thresholds["mid"])
select_example_years(waikato_FP,  thr_list[["waikato_FP"]]$thresholds["mid"])

select_example_years(napier_CD,   thr_list[["napier_CD"]]$thresholds["mid"])
select_example_years(napier_FP,   thr_list[["napier_FP"]]$thresholds["mid"])

select_example_years(northland_CD,thr_list[["northland_CD"]]$thresholds["mid"])
select_example_years(northland_FP,thr_list[["northland_FP"]]$thresholds["mid"])

select_example_years(milford_CD,  thr_list[["milford_CD"]]$thresholds["mid"])
select_example_years(milford_FP,  thr_list[["milford_FP"]]$thresholds["mid"])

select_example_years(waikato_CD,  thr_list[["waikato_CD"]]$thresholds["onein360"])
select_example_years(waikato_FP,  thr_list[["waikato_FP"]]$thresholds["onein360"])

select_example_years(napier_CD,   thr_list[["napier_CD"]]$thresholds["onein360"])
select_example_years(napier_FP,   thr_list[["napier_FP"]]$thresholds["onein360"])

select_example_years(northland_CD,thr_list[["northland_CD"]]$thresholds["onein360"])
select_example_years(northland_FP,thr_list[["northland_FP"]]$thresholds["onein360"])

select_example_years(milford_CD,  thr_list[["milford_CD"]]$thresholds["onein360"])
select_example_years(milford_FP,  thr_list[["milford_FP"]]$thresholds["onein360"])


day_cols <- as.character(1:360)

annual_rain <- function(df) {
  day_cols <- as.character(1:360)
  
  df %>%
    mutate(
      AnnualRain = rowSums(across(all_of(day_cols)), na.rm = TRUE)
    )
}

waikato_CD  <- annual_rain(waikato_CD)
waikato_FP  <- annual_rain(waikato_FP)
milford_CD  <- annual_rain(milford_CD)
milford_FP  <- annual_rain(milford_FP)
napier_CD   <- annual_rain(napier_CD)
napier_FP   <- annual_rain(napier_FP)
northland_CD<- annual_rain(northland_CD)
northland_FP<- annual_rain(northland_FP)

summary(waikato_CD$AnnualRain)

make_example_year_row <- function(df, exceedances, region, scenario, threshold, category, year) {
  
  if (!"AnnualRain" %in% colnames(df)) {
    day_cols <- as.character(1:360)
    df <- df %>%
      mutate(AnnualRain = rowSums(across(all_of(day_cols)), na.rm = TRUE))
  }
  
  n_years <- nrow(df)
  
  rx_stats <- df %>%
    arrange(desc(RX1day)) %>%
    mutate(
      RX1day_rank = row_number(),
      RX1day_percentile = (1 - (RX1day_rank - 1) / (n_years - 1)) * 100
    )
  
  yr_row <- df %>% filter(Year == year)
  yr_rx  <- rx_stats %>% filter(Year == year)
  
  data.frame(
    Region = region,
    Scenario = scenario,
    Threshold = threshold,
    `Example Type` = category,
    Year = year,
    `Exceedance Days` = exceedances[as.character(year)],
    `RX1day (mm)` = yr_row$RX1day,
    `RX1day Rank` = yr_rx$RX1day_rank,
    `RX1day Percentile` = yr_rx$RX1day_percentile,
    `Annual Rainfall (mm)` = yr_row$AnnualRain,
    stringsAsFactors = FALSE,
    check.names = FALSE)
}

build_example_year_table <- function(df, exceedances, region, scenario, threshold_label, threshold_value) {
  
  yrs <- select_example_years(df, threshold_value)
  
  bind_rows(
    make_example_year_row(df, exceedances, region, scenario, threshold_label,
                          "Muted", yrs$muted),
    
    make_example_year_row(df, exceedances, region, scenario, threshold_label,
                          "Single Exceedance", yrs$single),
    
    make_example_year_row(df, exceedances, region, scenario, threshold_label,
                          "High Exceedance", yrs$high)
  )
}

table_example_years_all_mid <- bind_rows(
  
  build_example_year_table(waikato_CD, waikato_cd_mid,
                           "Waikato", "Current Day", "Mid RX1day",
                           thr_list[["waikato_CD"]]$thresholds["mid"]),
  
  build_example_year_table(
    waikato_FP, waikato_fp_mid,
    "Waikato", "Future Projection", "Mid RX1day",
    thr_list[["waikato_FP"]]$thresholds["mid"]),
  
  build_example_year_table(napier_CD, napier_cd_mid,
                           "Napier", "Current Day", "Mid RX1day",
                           thr_list[["napier_CD"]]$thresholds["mid"]),
  
  build_example_year_table(napier_FP, napier_fp_mid,
                           "Napier", "Future Projection", "Mid RX1day",
                           thr_list[["napier_FP"]]$thresholds["mid"]),
  
  build_example_year_table(northland_CD, northland_cd_mid,
                           "Northland", "Current Day", "Mid RX1day",
                           thr_list[["northland_CD"]]$thresholds["mid"]),
  
  build_example_year_table(northland_FP, northland_fp_mid,
                           "Northland", "Future Projection", "Mid RX1day",
                           thr_list[["northland_FP"]]$thresholds["mid"]),
  
  build_example_year_table(milford_CD, milford_cd_mid,
                           "Milford Sound", "Current Day", "Mid RX1day",
                           thr_list[["milford_CD"]]$thresholds["mid"]),
  
  build_example_year_table(milford_FP, milford_fp_mid,
                           "Milford Sound", "Future Projection", "Mid RX1day",
                           thr_list[["milford_FP"]]$thresholds["mid"]))
table_example_years_all_mid


table_example_years_all_onein360 <- bind_rows(
  
  build_example_year_table(waikato_CD, waikato_cd_360,
                           "Waikato", "Current Day", "1 in 360",
                           thr_list[["waikato_CD"]]$thresholds["onein360"]),
  
  build_example_year_table(waikato_FP, waikato_fp_360,
                           "Waikato", "Future Projection", "1 in 360",
                           thr_list[["waikato_FP"]]$thresholds["onein360"]),
  
  build_example_year_table(napier_CD, napier_cd_360,
                           "Napier", "Current Day", "1 in 360",
                           thr_list[["napier_CD"]]$thresholds["onein360"]),
  
  build_example_year_table(napier_FP, napier_fp_360,
                           "Napier", "Future Projection", "1 in 360",
                           thr_list[["napier_FP"]]$thresholds["onein360"]),
  
  build_example_year_table(northland_CD, northland_cd_360,
                           "Northland", "Current Day", "1 in 360",
                           thr_list[["northland_CD"]]$thresholds["onein360"]),
  
  build_example_year_table(northland_FP, northland_fp_360,
                           "Northland", "Future Projection", "1 in 360",
                           thr_list[["northland_FP"]]$thresholds["onein360"]),
  
  build_example_year_table(milford_CD, milford_cd_360,
                           "Milford Sound", "Current Day", "1 in 360",
                           thr_list[["milford_CD"]]$thresholds["onein360"]),
  
  build_example_year_table(milford_FP, milford_fp_360,
                           "Milford Sound", "Future Projection", "1 in 360",
                           thr_list[["milford_FP"]]$thresholds["onein360"]))

table_example_years_all_onein360

table_example_years_appendix <- bind_rows(
  table_example_years_all_mid,
  table_example_years_all_onein360
) %>%
  arrange(Region, Scenario, Threshold, `Example Type`)

write.csv(
  table_example_years_appendix,
  file = file.path(model_data_dir, "example_years_table_all.csv"),
  row.names = FALSE
)

write.csv(
  table_example_years_all_onein360,
  file = file.path(model_data_dir, "example_years_table_onein360.csv"),
  row.names = FALSE
)

write.csv(
  table_example_years_all_mid,
  file = file.path(model_data_dir, "example_years_table_mid.csv"),
  row.names = FALSE
)


# Addressing Suggestions --------------------------------------------------

# Fixed-threshold histograms (current thresholds applied to Future) -------

build_hist_df_fixedCD_all <- function(region_name, thr_list, df_CD, df_FP) {
  
  thr_CD <- thr_list[[paste0(region_name, "_CD")]]$thresholds
  
  hist_df <- bind_rows(data.frame(days = count_exceedances_per_year(df_CD, thr_CD["mid"]),Threshold = "> Mid RX1day", Period = "Current Day"),
                       data.frame(days = count_exceedances_per_year(df_CD, thr_CD["onein360"]),Threshold = "> 1-in-360 threshold", Period = "Current Day"),
                       data.frame(days = count_exceedances_per_year(df_CD, thr_CD["mean"]), Threshold = "> Mean RX1day", Period = "Current Day"),
                       
                       # Future Projection (CD thresholds on FP data)
                       data.frame(days = count_exceedances_per_year(df_FP, thr_CD["mid"]), Threshold = "> Mid RX1day", Period = "Future Projection"),
                       data.frame(days = count_exceedances_per_year(df_FP, thr_CD["onein360"]), Threshold = "> 1-in-360 threshold", Period = "Future Projection"),
                       data.frame(days = count_exceedances_per_year(df_FP, thr_CD["mean"]), Threshold = "> Mean RX1day", Period = "Future Projection"))
  
  hist_df$Threshold <- factor(hist_df$Threshold, levels = c( "> Mid RX1day", "> 1-in-360 threshold", "> Mean RX1day"))
  
  hist_df_prop <- hist_df %>%
    count(Period, Threshold, days) %>%
    group_by(Period, Threshold) %>%
    mutate(prop_years = n / sum(n)) %>%
    ungroup()
  
  hist_df_prop$Region <- tools::toTitleCase(region_name)
  hist_df_prop}

regions_mod <- c("waikato", "napier", "northland", "milford")

hist_fixedCD_all <- lapply(regions_mod, function(reg) {build_hist_df_fixedCD_all(region_name = reg, thr_list = thr_list, df_CD = get(paste0(reg, "_CD")), df_FP = get(paste0(reg, "_FP")))})

names(hist_fixedCD_all) <- regions_mod
hist_plots_fixedCD <- lapply(hist_fixedCD_all, plot_hist_exceedances)

hist_plots_fixedCD$waikato
hist_plots_fixedCD$milford

for (reg in regions_mod) {save_plot(hist_plots_fixedCD[[reg]], paste0(reg, "_exceedance_fixedCD_allthresholds_histogram.png"), width = fig_width_standard, height = fig_height_standard)}

#checking it is working 
hist_original <- build_hist_df("waikato", thr_list, waikato_CD, waikato_FP)
hist_fixed    <- build_hist_df_fixedCD_all("waikato", thr_list, waikato_CD, waikato_FP)

# Compare CD rows only
hist_original %>% filter(Period == "Current Day") %>% arrange(Threshold, days)
hist_fixed %>% filter(Period == "Current Day") %>%arrange(Threshold, days)
#The Current day proportions are the same as expected. 

# Summary table for histogram  --------------------------------------------

#Change in mean number of exceedances
build_fixedCD_summary_table <- function(hist_df_prop) {
  hist_df_prop %>%
    group_by(Region, Threshold, Period) %>%
    summarise(Mean_exceedances = weighted.mean(days, prop_years), .groups = "drop")}

fixedCD_changeinmean_table <- fixedCD_summary_all %>%
  pivot_wider(names_from  = Period, values_from = Mean_exceedances) %>%
  rename(`Mean Exceedance CD` = `Current Day`, `Mean Exceedance FP` = `Future Projection`) %>%
  mutate(Mean_change = `Mean Exceedance FP` - `Mean Exceedance CD`)

fixedCD_changeinmean_table

#Cumulative proportions. 
calc_cumulative_props <- function(exceed_vec) {
  max_k <- max(exceed_vec, na.rm = TRUE)
  n_years <- length(exceed_vec)
  data.frame(k = 0:max_k, prop_ge_k = sapply(0:max_k, function(k) {sum(exceed_vec >= k, na.rm = TRUE) / n_years}))}

build_cumulative_exceed_table <- function(region_name, threshold_name, thr_list, df_CD, df_FP) {
  thr_value <- thr_list[[paste0(region_name, "_CD")]]$thresholds[threshold_name]
  exc_CD <- count_exceedances_per_year(df_CD, thr_value)
  exc_FP <- count_exceedances_per_year(df_FP, thr_value)
  max_k <- max(c(exc_CD, exc_FP), na.rm = TRUE)
  k_seq <- 0:max_k
  cum_CD <- data.frame(k = k_seq, prop_CD = sapply(k_seq, function(k) mean(exc_CD >= k)))
  cum_FP <- data.frame(k = k_seq, prop_FP = sapply(k_seq, function(k) mean(exc_FP >= k)))
  summary_df <- left_join(cum_CD, cum_FP, by = "k") %>%
    mutate(prop_CD = replace_na(prop_CD, 0), prop_FP = replace_na(prop_FP, 0),
           FP_to_CD = ifelse(prop_CD == 0, NA, prop_FP / prop_CD),
           Region = tools::toTitleCase(region_name),
           Threshold = threshold_name)
  summary_df}

regions_mod <- c("waikato", "napier", "northland", "milford")
thresholds <- c("mid", "mean", "onein360")

cumulative_proportion_table <- bind_rows(lapply(regions_mod, function(reg) {
  bind_rows(lapply(thresholds, function(thr) {
    build_cumulative_exceed_table(
      region_name = reg,
      threshold_name = thr,
      thr_list = thr_list,
      df_CD = get(paste0(reg, "_CD")),
      df_FP = get(paste0(reg, "_FP")))}))
}))

threshold_labels_cumulative_proportion_table <- c(mid = "> Mid RX1day", mean = "> Mean RX1day", onein360 = "> 1-in-360")

cumulative_proportion_table <- cumulative_proportion_table %>%
  mutate(Threshold = threshold_labels_cumulative_proportion_table[Threshold],
         prop_CD  = round(prop_CD, 4), prop_FP  = round(prop_FP, 4), FP_to_CD = round(FP_to_CD, 4)) %>%
  rename(`Proportion CD` = prop_CD, `Proportion FP` = prop_FP,`FP / CD` = FP_to_CD, `Exceedances ≥ k` = k) %>%
  select(Region, Threshold,`Exceedances ≥ k`,`Proportion CD`, `Proportion FP`,`FP / CD`)

cumulative_proportion_table

# Save cumulative proportion table as CSV
write.csv(
  cumulative_proportion_table,
  file = file.path(model_data_dir, "cumulative_proportion_table.csv"),
  row.names = FALSE
)

#Updating boxplots with fixed thresholds from current day. 
plot_rx1day_vs_exceedance_fixedCD <- function(data, exceedances, title_text, baseline_threshold = NULL, exceed_max = NULL) {
  
  if (!is.null(baseline_threshold)) {
    exceedances <- count_exceedances_per_year(data, baseline_threshold)
  }
  
  df_plot <- data %>%
    mutate(exceedances = exceedances) %>%
    select(Year, RX1day, exceedances)
  
  # Compute max RX1day for plotting
  rx1day_max <- max(df_plot$RX1day, na.rm = TRUE) * 1.1
  
  # If exceed_max is not given, compute from data
  if (is.null(exceed_max)) {
    exceed_max <- max(df_plot$exceedances, na.rm = TRUE)
  }
  
  # Calculate percent of points per x-axis bin
  percent_per_bin <- df_plot %>%
    group_by(exceedances) %>%
    summarise(percent = n() / nrow(df_plot) * 100, .groups = "drop")
  
  ggplot(df_plot, aes(x = factor(exceedances), y = RX1day, fill = factor(exceedances))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(data = function(d) {
      box_stats <- boxplot.stats(d$RX1day)$stats
      d %>% filter(RX1day < box_stats[1] | RX1day > box_stats[5])
    }, width = 0.15, alpha = 0.4) +
    geom_text(
      data = percent_per_bin,
      aes(
        x = factor(exceedances),
        y = rx1day_max * 1.02,
        label = ifelse(percent < 0.1, "<0.1%", paste0(round(percent, 1), "%"))
      ),
      inherit.aes = FALSE,
      vjust = 0,
      size = 2.3,
      hjust = 0.5
    ) +
    theme_thesis +
    labs(x = "Number of exceedances per year", y = "RX1day (mm)", title = title_text) +
    guides(fill = "none") +
    scale_y_continuous(breaks = seq(0, rx1day_max, by = 100)) +
    coord_cartesian(ylim = c(0, rx1day_max * 1.1), clip = "off") +
    scale_x_discrete(limits = as.character(0:exceed_max))
}
region_exceed_max <- list()

for (reg in regions_mod) {
  all_exceedances <- c(
    get(paste0(reg, "_cd_mid")),
    get(paste0(reg, "_cd_360")),
    get(paste0(reg, "_fp_mid")),
    get(paste0(reg, "_fp_360"))
  )
  
  region_exceed_max[[reg]] <- max(all_exceedances, na.rm = TRUE)
}


#checking it is working
boxplot_fp_mid_waikato_fixedCD <- plot_rx1day_vs_exceedance_fixedCD(
  waikato_FP, waikato_fp_mid, 
  "Mid RX1day - Future Projection (Current Threshold)",
  baseline_threshold = thr_list[["waikato_CD"]]$thresholds["mid"])

boxplot_fp_mid_waikato_fixedCD

for (reg in regions_mod) {
  box_cd_mid_fixedCD  <- plot_rx1day_vs_exceedance_fixedCD(
    get(paste0(reg, "_CD")), 
    get(paste0(reg, "_cd_mid")), 
    "Mid RX1day - Current Day",
    baseline_threshold = thr_list[[paste0(reg, "_CD")]]$thresholds["mid"],
    exceed_max = region_exceed_max[[reg]]
  )
  
  box_fp_mid_fixedCD  <- plot_rx1day_vs_exceedance_fixedCD(
    get(paste0(reg, "_FP")), 
    get(paste0(reg, "_fp_mid")), 
    "Mid RX1day - Future Projection",
    baseline_threshold = thr_list[[paste0(reg, "_CD")]]$thresholds["mid"],
    exceed_max = region_exceed_max[[reg]]
  )
  
  box_cd_360_fixedCD  <- plot_rx1day_vs_exceedance_fixedCD(
    get(paste0(reg, "_CD")), 
    get(paste0(reg, "_cd_360")), 
    "1-in-360 - Current Day",
    baseline_threshold = thr_list[[paste0(reg, "_CD")]]$thresholds["onein360"],
    exceed_max = region_exceed_max[[reg]]
  )
  
  box_fp_360_fixedCD  <- plot_rx1day_vs_exceedance_fixedCD(
    get(paste0(reg, "_FP")), 
    get(paste0(reg, "_fp_360")), 
    "1-in-360 - Future Projection",
    baseline_threshold = thr_list[[paste0(reg, "_CD")]]$thresholds["onein360"],
    exceed_max = region_exceed_max[[reg]]
  )
  
  # Combine into 2x2 layout
  combined_boxplot_fixedCD <- (box_cd_mid_fixedCD | box_fp_mid_fixedCD) /
    (box_cd_360_fixedCD | box_fp_360_fixedCD) +
    plot_annotation(title = region_labels[reg],
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold")))
  
  # Save the combined plot
  save_plot(combined_boxplot_fixedCD, paste0(reg, "_rx1day_all_boxplots_fixedCD.png"), width = fig_width_hoz_half, height = fig_height_hoz_half)
}
