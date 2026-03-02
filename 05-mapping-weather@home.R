# -------------------------------------------------------------------------
# 05-mapping-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Mapping the probability ratio (future/current) for years with
# greater than or equal to 4 exceedance days at each weather@home grid cell.
# -------------------------------------------------------------------------

# Packages -----------------------------------------------------------------
# install.packages(c("ggplot2", "viridis"))
library(ggplot2)
library(viridis)

# Input/output paths --------------------------------------------------------
weatherathome_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"

input_file <- file.path(
  weatherathome_dir,
  "weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv"
)

output_png <- file.path(
  weatherathome_dir,
  "weather@home_probability_ratio_ge4_map.png"
)

output_csv <- file.path(
  weatherathome_dir,
  "weather@home_probability_ratio_ge4_map_data.csv"
)

# Read outputs from Script 04 ----------------------------------------------
if (!file.exists(input_file)) {
  stop("Input file not found. Run 04-processing-weather@home.R first.")
}

grid_results <- read.csv(input_file)

required_cols <- c(
  "global_longitude0",
  "global_latitude0",
  "probability_ratio_ge4_future_over_current"
)

missing_cols <- setdiff(required_cols, names(grid_results))
if (length(missing_cols) > 0) {
  stop(
    sprintf(
      "Input file is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Keep finite values for colour scaling; track Inf/NA for diagnostics.
map_data <- data.frame(
  lon = grid_results$global_longitude0,
  lat = grid_results$global_latitude0,
  probability_ratio_ge4 = grid_results$probability_ratio_ge4_future_over_current
)

map_data$ratio_class <- ifelse(
  is.na(map_data$probability_ratio_ge4),
  "NA (current=0, future=0)",
  ifelse(
    is.infinite(map_data$probability_ratio_ge4),
    "Inf (current=0, future>0)",
    "Finite"
  )
)

finite_map_data <- map_data[is.finite(map_data$probability_ratio_ge4), ]

if (nrow(finite_map_data) == 0) {
  stop("No finite probability-ratio values available to plot.")
}

# Optional clipping to make map colours easier to interpret.
# This keeps the plotted scale stable if a few cells are very large.
q99 <- as.numeric(quantile(finite_map_data$probability_ratio_ge4, probs = 0.99, na.rm = TRUE))
finite_map_data$probability_ratio_ge4_plot <- pmin(finite_map_data$probability_ratio_ge4, q99)

# Save plotting data for reproducibility
write.csv(map_data, output_csv, row.names = FALSE)

# Plot ---------------------------------------------------------------------
p_ge4_ratio <- ggplot(finite_map_data, aes(x = lon, y = lat, fill = probability_ratio_ge4_plot)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_viridis(
    option = "magma",
    name = "Probability\nratio",
    direction = 1
  ) +
  labs(
    title = "weather@home: Probability ratio for >=4 exceedance days",
    subtitle = "Future (3k warmer) / Current decade, by model grid cell",
    x = "Longitude",
    y = "Latitude",
    caption = paste0(
      "Finite values are plotted. Values above the 99th percentile are clipped at ",
      round(q99, 2),
      "."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(
  filename = output_png,
  plot = p_ge4_ratio,
  width = 8,
  height = 7,
  dpi = 300
)

# Quick diagnostics ---------------------------------------------------------
finite_count <- sum(is.finite(map_data$probability_ratio_ge4))
inf_count <- sum(is.infinite(map_data$probability_ratio_ge4))
na_count <- sum(is.na(map_data$probability_ratio_ge4))

cat("Finite cell count:", finite_count, "\n")
cat("Infinite cell count:", inf_count, "\n")
cat("NA cell count:", na_count, "\n")
cat("Finite ratio min:", round(min(finite_map_data$probability_ratio_ge4, na.rm = TRUE), 4), "\n")
cat("Finite ratio median:", round(median(finite_map_data$probability_ratio_ge4, na.rm = TRUE), 4), "\n")
cat("Finite ratio max:", round(max(finite_map_data$probability_ratio_ge4, na.rm = TRUE), 4), "\n")
cat("Saved map to:", output_png, "\n")
