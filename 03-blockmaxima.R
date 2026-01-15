# -------------------------------------------------------------------------
# 03-blockmaxima.R
# -------------------------------------------------------------------------
# Nov 2025
# comparing blockmaxima to all in hydrological year.
# -------------------------------------------------------------------------

# Filter Tairua Forest ----------------------------------------------------
tairua_forest_df <- combined_df %>%
  filter(station == "Tairua Forest") %>%
  select(observation_date, rainfall_mm, hydro_year)

# RX1day per hydro year ---------------------------------------------------
tairua_forest_rx1day_df <- tairua_forest_df %>%
  group_by(hydro_year) %>%
  summarise(
    RX1day = max(rainfall_mm, na.rm = TRUE),
    .groups = "drop"
  )

# threshold values --------------------------------------------------------
tairua_forest_minRX1day <- min(tairua_forest_rx1day_df$RX1day, na.rm = TRUE)
tairua_forest_meanRX1day <- mean(tairua_forest_rx1day_df$RX1day, na.rm = TRUE)
tairua_forest_midRX1day <- (tairua_forest_minRX1day + tairua_forest_meanRX1day) / 2
tairua_forest_maxRX1day <- max(tairua_forest_rx1day_df$RX1day, na.rm = TRUE)

# Proportion of days above threshold — Tairua Forest ---------------------
# selecting only daily rainfall column
tairua_forest_daily_vec <- tairua_forest_df$rainfall_mm

# number of daily exceedances for each threshold
tairua_forest_n_daily_exceed_min  <- sum(tairua_forest_daily_vec > tairua_forest_minRX1day,  na.rm = TRUE)
tairua_forest_n_daily_exceed_mean <- sum(tairua_forest_daily_vec > tairua_forest_meanRX1day, na.rm = TRUE)
tairua_forest_n_daily_exceed_mid  <- sum(tairua_forest_daily_vec > tairua_forest_midRX1day,  na.rm = TRUE)
tairua_forest_n_daily_exceed_max  <- sum(tairua_forest_daily_vec > tairua_forest_maxRX1day,  na.rm = TRUE)

# total number of valid daily observations
tairua_forest_n_days <- sum(!is.na(tairua_forest_daily_vec))

# proportion of days above each threshold
tairua_forest_prop_exceed_min  <- tairua_forest_n_daily_exceed_min  / tairua_forest_n_days
tairua_forest_prop_exceed_mean <- tairua_forest_n_daily_exceed_mean / tairua_forest_n_days
tairua_forest_prop_exceed_mid  <- tairua_forest_n_daily_exceed_mid  / tairua_forest_n_days
tairua_forest_prop_exceed_max  <- tairua_forest_n_daily_exceed_max  / tairua_forest_n_days

print(tairua_forest_prop_exceed_min)
print(tairua_forest_prop_exceed_mid)
print(tairua_forest_prop_exceed_mean)
print(tairua_forest_prop_exceed_max)

# 1-in-360-day threshold
tairua_forest_prop_exceed_onein360 <- 1 / 360
tairua_forest_onein360_quantile <- 1 - tairua_forest_prop_exceed_onein360

tairua_forest_threshold_onein360 <- quantile(
  tairua_forest_daily_vec,
  probs = tairua_forest_onein360_quantile,
  na.rm = TRUE,
  type = 7
)
print(tairua_forest_threshold_onein360)

# labels for plots including threshold value + percent exceedance
labels_tairua_forest <- c(
  "Min RX1day" = paste0(
    "Min RX1day (",
    round(tairua_forest_minRX1day, 1), " mm, ",
    round(tairua_forest_prop_exceed_min * 100, 2), "%)"
  ),
  "Mid RX1day" = paste0(
    "Mid RX1day (",
    round(tairua_forest_midRX1day, 1), " mm, ",
    round(tairua_forest_prop_exceed_mid * 100, 2), "%)"
  ),
  "Mean RX1day" = paste0(
    "Mean RX1day (",
    round(tairua_forest_meanRX1day, 1), " mm, ",
    round(tairua_forest_prop_exceed_mean * 100, 2), "%)"
  ),
  "1 in 360 day threshold" = paste0(
    "1 in 360 day (",
    round(tairua_forest_threshold_onein360, 1), " mm, ",
    round(tairua_forest_prop_exceed_onein360 * 100, 2), "%)"
  ),
  "Max RX1day" = paste0(
    "Max RX1day (",
    round(tairua_forest_maxRX1day, 1), " mm, ",
    round(tairua_forest_prop_exceed_max * 100, 2), "%)"
  )
)

labels_tairua_forest

# RX1day time series plot -------------------------------------------------
p_tairua_forest <- ggplot(tairua_forest_rx1day_df, aes(x = hydro_year, y = RX1day)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  geom_hline(aes(yintercept = tairua_forest_minRX1day, colour = "Min RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = tairua_forest_meanRX1day, colour = "Mean RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = tairua_forest_midRX1day, colour = "Mid RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = tairua_forest_maxRX1day, colour = "Max RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = tairua_forest_threshold_onein360, colour = "1 in 360 day threshold"),
             linetype = "solid", size = 1.1) +
  scale_colour_manual(
    breaks = c(
      "Max RX1day",
      "Mean RX1day",
      "1 in 360 day threshold",
      "Mid RX1day",
      "Min RX1day"
    ),
    labels = labels_tairua_forest,
    values = c(
      "Min RX1day"             = "#56B4E9",  
      "Mid RX1day"             = "#0072B2", 
      "1 in 360 day threshold" = "#E69F00",  
      "Mean RX1day"            = "#009E73",  
      "Max RX1day"             = "#D55E00"   
    )
  ) +
  labs(
    title = "RX1day Time Series — Tairua Forest",
    x = "Hydrological Year (July–June)",
    y = "Maximum Daily Rainfall (mm)",
    colour = "Threshold"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_tairua_forest

# Filter Chiltern ----------------------------------------------------
chiltern_df <- combined_df %>%
  filter(station == "Chiltern") %>%
  select(observation_date, rainfall_mm, hydro_year)

# RX1day per hydro year ---------------------------------------------------
chiltern_rx1day_df <- chiltern_df %>%
  group_by(hydro_year) %>%
  summarise(
    RX1day = max(rainfall_mm, na.rm = TRUE),
    .groups = "drop"
  )

# threshold values --------------------------------------------------------
chiltern_minRX1day <- min(chiltern_rx1day_df$RX1day, na.rm = TRUE)
chiltern_meanRX1day <- mean(chiltern_rx1day_df$RX1day, na.rm = TRUE)
chiltern_midRX1day <- (chiltern_minRX1day + chiltern_meanRX1day) / 2
chiltern_maxRX1day <- max(chiltern_rx1day_df$RX1day, na.rm = TRUE)

# Proportion of days above threshold — Chiltern ---------------------
# selecting only daily rainfall column
chiltern_daily_vec <- chiltern_df$rainfall_mm

# number of daily exceedances for each threshold
chiltern_n_daily_exceed_min  <- sum(chiltern_daily_vec > chiltern_minRX1day,  na.rm = TRUE)
chiltern_n_daily_exceed_mean <- sum(chiltern_daily_vec > chiltern_meanRX1day, na.rm = TRUE)
chiltern_n_daily_exceed_mid  <- sum(chiltern_daily_vec > chiltern_midRX1day,  na.rm = TRUE)
chiltern_n_daily_exceed_max  <- sum(chiltern_daily_vec > chiltern_maxRX1day,  na.rm = TRUE)

# total number of valid daily observations
chiltern_n_days <- sum(!is.na(chiltern_daily_vec))

# proportion of days above each threshold
chiltern_prop_exceed_min  <- chiltern_n_daily_exceed_min  / chiltern_n_days
chiltern_prop_exceed_mean <- chiltern_n_daily_exceed_mean / chiltern_n_days
chiltern_prop_exceed_mid  <- chiltern_n_daily_exceed_mid  / chiltern_n_days
chiltern_prop_exceed_max  <- chiltern_n_daily_exceed_max  / chiltern_n_days

print(chiltern_prop_exceed_min)
print(chiltern_prop_exceed_mid)
print(chiltern_prop_exceed_mean)
print(chiltern_prop_exceed_max)

# 1-in-360-day threshold
chiltern_prop_exceed_onein360 <- 1 / 360
chiltern_onein360_quantile <- 1 - chiltern_prop_exceed_onein360

chiltern_threshold_onein360 <- quantile(
  chiltern_daily_vec,
  probs = chiltern_onein360_quantile,
  na.rm = TRUE,
  type = 7
)
print(chiltern_threshold_onein360)

# labels for plots including threshold value + percent exceedance
labels_chiltern <- c(
  "Min RX1day" = paste0(
    "Min RX1day (",
    round(chiltern_minRX1day, 1), " mm, ",
    round(chiltern_prop_exceed_min * 100, 2), "%)"
  ),
  "Mid RX1day" = paste0(
    "Mid RX1day (",
    round(chiltern_midRX1day, 1), " mm, ",
    round(chiltern_prop_exceed_mid * 100, 2), "%)"
  ),
  "Mean RX1day" = paste0(
    "Mean RX1day (",
    round(chiltern_meanRX1day, 1), " mm, ",
    round(chiltern_prop_exceed_mean * 100, 2), "%)"
  ),
  "1 in 360 day threshold" = paste0(
    "1 in 360 day (",
    round(chiltern_threshold_onein360, 1), " mm, ",
    round(chiltern_prop_exceed_onein360 * 100, 2), "%)"
  ),
  "Max RX1day" = paste0(
    "Max RX1day (",
    round(chiltern_maxRX1day, 1), " mm, ",
    round(chiltern_prop_exceed_max * 100, 2), "%)"
  )
)

labels_chiltern

# RX1day time series plot -------------------------------------------------
p_chiltern <- ggplot(chiltern_rx1day_df, aes(x = hydro_year, y = RX1day)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  geom_hline(aes(yintercept = chiltern_minRX1day, colour = "Min RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = chiltern_meanRX1day, colour = "Mean RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = chiltern_midRX1day, colour = "Mid RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = chiltern_maxRX1day, colour = "Max RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = chiltern_threshold_onein360, colour = "1 in 360 day threshold"),
             linetype = "solid", size = 1.1) +
  scale_colour_manual(
    breaks = c(
      "Max RX1day",
      "Mean RX1day",
      "1 in 360 day threshold",
      "Mid RX1day",
      "Min RX1day"
    ),
    labels = labels_chiltern,
    values = c(
      "Min RX1day"             = "#56B4E9",  
      "Mid RX1day"             = "#0072B2", 
      "1 in 360 day threshold" = "#E69F00",  
      "Mean RX1day"            = "#009E73",  
      "Max RX1day"             = "#D55E00"   
    )
  ) +
  labs(
    title = "RX1day Time Series — Chiltern",
    x = "Hydrological Year (July–June)",
    y = "Maximum Daily Rainfall (mm)",
    colour = "Threshold"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p_chiltern

# histogram ---------------------------------------------------------------

# Helper function 
count_exceedances_per_hydro_year <- function(df, threshold) {
  df %>%
    group_by(hydro_year) %>%
    summarise(
      days = sum(rainfall_mm > threshold, na.rm = TRUE),
      .groups = "drop"
    )
}

# Tairua Forest exceedance counts
tairua_mid_df <- count_exceedances_per_hydro_year(
  tairua_forest_df,
  tairua_forest_midRX1day
)

tairua_360_df <- count_exceedances_per_hydro_year(
  tairua_forest_df,
  tairua_forest_threshold_onein360
)

tairua_mean_df <- count_exceedances_per_hydro_year(
  tairua_forest_df,
  tairua_forest_meanRX1day
)

# Chiltern exceedance counts
chiltern_mid_df <- count_exceedances_per_hydro_year(
  chiltern_df,
  chiltern_midRX1day
)

chiltern_360_df <- count_exceedances_per_hydro_year(
  chiltern_df,
  chiltern_threshold_onein360
)

chiltern_mean_df <- count_exceedances_per_hydro_year(
  chiltern_df,
  chiltern_meanRX1day
)

# setting plot limits
x_lim <- c(-0.5, 12.5)
y_lim <- c(0, 0.6)

# Tauira Plots
p_tairua_mid <- ggplot(tairua_mid_df, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#0072B2"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> Mid RX1day",
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

p_tairua_360 <- ggplot(tairua_360_df, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#E69F00"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> 1-in-360 day threshold",
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

p_tairua_mean <- ggplot(tairua_mean_df, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#009E73"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> Mean RX1day",
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

# Combine — Tairua Forest
p_tairua_hist <- p_tairua_mid | p_tairua_360 | p_tairua_mean
p_tairua_hist


p_tairua_hist_final <- p_tairua_hist +
  plot_annotation(
    title = "Tairua Forest",
    theme = theme(
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      )
    )
  )
p_tairua_hist_final

# Chiltern 
p_chiltern_mid <- ggplot(chiltern_mid_df, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#0072B2"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> Mid RX1day",
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

p_chiltern_360 <- ggplot(chiltern_360_df, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#E69F00"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> 1-in-360 day threshold",
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

p_chiltern_mean <- ggplot(chiltern_mean_df, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#009E73"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> Mean RX1day",
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

# Combine — Chiltern
p_chiltern_hist <- p_chiltern_mid | p_chiltern_360 | p_chiltern_mean
p_chiltern_hist

p_chiltern_hist_final <- p_chiltern_hist +
  plot_annotation(
    title = "Chiltern",
    theme = theme(
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      )
    )
  )

p_chiltern_hist_final

