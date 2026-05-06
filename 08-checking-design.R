# -------------------------------------------------------------------------
# 08-checking-design.R
# -------------------------------------------------------------------------
# April 2026
# Aligning Extreme threshold with % of at least 4 excedance day in current day. 
# -------------------------------------------------------------------------

# Matching Rx1day Extreme threshold with 1.1% -----------------------------

rx1day_threshold_98.9_current <- apply(
  current_rx_array,
  c(1, 2),
  quantile,
  probs = 0.989,
  na.rm = TRUE,
  type = 7
)

rx1day_threshold_98.9_future <- apply(
  future_rx_array,
  c(1, 2),
  quantile,
  probs = 0.989,
  na.rm = TRUE,
  type = 7
)

current_prop_rx1day_98.9 <- calc_rx1day_top10_proportion(
  current_rx_array,
  rx1day_threshold_98.9_current
)
future_prop_rx1day_98.9 <- calc_rx1day_top10_proportion(
  future_rx_array,
  rx1day_threshold_98.9_future
)

probability_ratio_rx1day_98.9 <- calc_probability_ratio(
  current_prop_rx1day_98.9,
  future_prop_rx1day_98.9
)

dim(current_prop_rx1day_98.9)
dim(future_prop_rx1day_98.9)
# both should be 44 X 44


# Matching Rx1day Extreme threshold with 3.4% -----------------------------

rx1day_threshold_96.6_current <- apply(
  current_rx_array,
  c(1, 2),
  quantile,
  probs = 0.966,
  na.rm = TRUE,
  type = 7
)

rx1day_threshold_96.6_future <- apply(
  future_rx_array,
  c(1, 2),
  quantile,
  probs = 0.966,
  na.rm = TRUE,
  type = 7
)

current_prop_rx1day_96.6 <- calc_rx1day_top10_proportion(
  current_rx_array,
  rx1day_threshold_96.6_current
)
future_prop_rx1day_96.6 <- calc_rx1day_top10_proportion(
  future_rx_array,
  rx1day_threshold_96.6_future
)

probability_ratio_rx1day_96.6 <- calc_probability_ratio(
  current_prop_rx1day_96.6,
  future_prop_rx1day_96.6
)

dim(current_prop_rx1day_96.6)
dim(future_prop_rx1day_96.6)
# both should be 44 X 44
