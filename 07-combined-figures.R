# -------------------------------------------------------------------------
# 07-combined-figures.R
# -------------------------------------------------------------------------
# Mar 2026
# Find the coordinates (row/column and optional lat/lon labels) for a
# specific weather@home grid cell using the land-sea mask value, then pull
# the matching value from the joint probability-ratio grid.
# -------------------------------------------------------------------------

# Setting input paths and target value ------------------------------------
land_mask_file <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/Land-Sea Mask for Weather@home Data.csv"
ratio_grid_file <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv"

target_value <- 294.9545
value_tolerance <- 1e-6

# Helper functions ---------------------------------------------------------
read_grid_csv <- function(file_path) {
  raw <- utils::read.csv(
    file = file_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  first_col <- raw[[1]]
  first_col_char <- trimws(as.character(first_col))
  suppressWarnings(first_col_num <- as.numeric(first_col_char))
  
  is_simple_index <- all(!is.na(first_col_num)) &&
    length(first_col_num) == nrow(raw) &&
    identical(first_col_num, seq_len(nrow(raw)))
  
  has_missing_or_blank <- any(
    is.na(first_col) |
      first_col_char == "" |
      tolower(first_col_char) %in% c("na", "nan")
  )
  
  has_unique_labels <- !anyDuplicated(first_col_char)
  
  use_first_col_as_rownames <- !is_simple_index && !has_missing_or_blank && has_unique_labels
  
  if (use_first_col_as_rownames) {
    rownames(raw) <- first_col_char
    raw <- raw[, -1, drop = FALSE]
  } else if (is_simple_index) {
    raw <- raw[, -1, drop = FALSE]
  }
  
  numeric_raw <- suppressWarnings(as.data.frame(lapply(raw, as.numeric), check.names = FALSE))
  grid_matrix <- as.matrix(numeric_raw)
  
  if (!is.null(rownames(raw))) {
    rownames(grid_matrix) <- rownames(raw)
  }
  colnames(grid_matrix) <- colnames(raw)
  
  grid_matrix
}

parse_numeric_labels <- function(labels) {
  if (is.null(labels)) {
    return(numeric(0))
  }
  suppressWarnings(as.numeric(labels))
}

find_matching_cells <- function(grid_matrix, value, tolerance = 1e-6) {
  which(abs(grid_matrix - value) <= tolerance, arr.ind = TRUE)
}

# Loading grid files -------------------------------------------------------
land_mask_matrix <- read_grid_csv(land_mask_file)
ratio_grid_matrix <- read_grid_csv(ratio_grid_file)

# Checking dimensions ------------------------------------------------------
if (!all(dim(land_mask_matrix) == dim(ratio_grid_matrix))) {
  stop(
    "Grid dimensions do not match.\n",
    "Land mask dimensions: ", paste(dim(land_mask_matrix), collapse = " x "), "\n",
    "Ratio grid dimensions: ", paste(dim(ratio_grid_matrix), collapse = " x "), "\n",
    "Check whether one file needs transposing or different row/column handling."
  )
}

# Finding matching grid cell(s) -------------------------------------------
match_indices <- find_matching_cells(land_mask_matrix, target_value, value_tolerance)

if (nrow(match_indices) == 0) {
  stop(
    "No grid cell found in the land mask for value ",
    target_value,
    " within tolerance ",
    value_tolerance,
    "."
  )
}

# Building output table ----------------------------------------------------
lat_labels <- rownames(land_mask_matrix)
lon_labels <- colnames(land_mask_matrix)

lat_values <- parse_numeric_labels(lat_labels)
lon_values <- parse_numeric_labels(lon_labels)

matched_cells <- data.frame(
  hit_number = seq_len(nrow(match_indices)),
  row_index = match_indices[, "row"],
  col_index = match_indices[, "col"],
  land_mask_value = land_mask_matrix[match_indices],
  ratio_grid_value = ratio_grid_matrix[match_indices],
  latitude = if (length(lat_values) > 0) lat_values[match_indices[, "row"]] else NA_real_,
  longitude = if (length(lon_values) > 0) lon_values[match_indices[, "col"]] else NA_real_
)

print(matched_cells)

# Saving output ------------------------------------------------------------
output_file <- "matched_grid_cells_for_294_9545.csv"
utils::write.csv(matched_cells, output_file, row.names = FALSE)

message("Saved matched grid cells to: ", output_file)

if (all(is.na(matched_cells$latitude)) || all(is.na(matched_cells$longitude))) {
  message(
    "Latitude/longitude labels were not detected from row/column names. ",
    "Use row_index and col_index from the output file as the exact cell coordinates."
  )
}
