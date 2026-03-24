# 07-combined-figures.R
#
# Purpose:
# Given a known value in the land-sea mask (e.g., 294.9545), identify the
# corresponding grid-cell row/column and extract that same cell from a second
# gridded CSV file. If latitude/longitude labels are present in row/column names,
# those coordinates are also returned.

options(stringsAsFactors = FALSE)

# ----------------------------- User inputs -----------------------------------
land_mask_path <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/Land-Sea Mask for Weather@home Data.csv"
ratio_grid_path <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv"

# The value you already know appears in the land mask
target_value <- 294.9545

# Tolerance for floating-point comparisons
value_tolerance <- 1e-6
# -----------------------------------------------------------------------------

read_grid_csv <- function(path) {
  # Read with names preserved (important when lon/lat are column headers)
  raw <- read.csv(path, check.names = FALSE)

  # Heuristic: if first column looks like a row-label/index column (often lat),
  # move it to row names before numeric conversion.
  first_col <- raw[[1]]
  first_col_name <- names(raw)[1]

  use_first_col_as_rownames <- FALSE

  # If the first column is non-numeric, it's very likely a label column.
  if (!is.numeric(first_col)) {
    suppressWarnings(first_col_num <- as.numeric(first_col))
    if (any(!is.na(first_col_num))) {
      # Mixed numeric/character: still likely labels
      use_first_col_as_rownames <- TRUE
    } else {
      use_first_col_as_rownames <- TRUE
    }
  } else {
    # If numeric but not a simple 1..n sequence, it could be latitude labels.
    seq_like <- all(!is.na(first_col)) && identical(first_col, seq_len(nrow(raw)))
    if (!seq_like) use_first_col_as_rownames <- TRUE
  }

  if (use_first_col_as_rownames) {
    rownames(raw) <- as.character(first_col)
    raw <- raw[, -1, drop = FALSE]
  }

  # Convert all remaining values to numeric.
  num <- suppressWarnings(as.data.frame(lapply(raw, as.numeric), check.names = FALSE))
  mat <- as.matrix(num)

  # Preserve row/column names where available.
  rownames(mat) <- if (!is.null(rownames(raw))) rownames(raw) else rownames(mat)
  colnames(mat) <- colnames(raw)

  list(
    matrix = mat,
    original_first_col_name = first_col_name,
    used_first_col_as_rownames = use_first_col_as_rownames
  )
}

parse_numeric_labels <- function(x) {
  if (is.null(x)) return(rep(NA_real_, 0))
  suppressWarnings(as.numeric(x))
}

find_grid_cells_by_value <- function(mat, value, tol = 1e-6) {
  idx <- which(abs(mat - value) <= tol, arr.ind = TRUE)
  if (nrow(idx) == 0) return(idx)
  idx
}

# Read both grids
land <- read_grid_csv(land_mask_path)
ratio <- read_grid_csv(ratio_grid_path)
land_mat <- land$matrix
ratio_mat <- ratio$matrix

# Basic shape check
if (!all(dim(land_mat) == dim(ratio_mat))) {
  stop(
    "Grid dimensions do not match.\n",
    "Land mask dimensions: ", paste(dim(land_mat), collapse = " x "), "\n",
    "Ratio grid dimensions: ", paste(dim(ratio_mat), collapse = " x "), "\n",
    "You may need to align/transpose one dataset before matching cells."
  )
}

# Find all matching cells in land mask
hits <- find_grid_cells_by_value(land_mat, target_value, tol = value_tolerance)

if (nrow(hits) == 0) {
  stop(
    "No cells found in land mask with value ", target_value,
    " within tolerance ", value_tolerance, "."
  )
}

# Pull coordinates (if encoded in names)
lat_labels <- rownames(land_mat)
lon_labels <- colnames(land_mat)
lat_values <- parse_numeric_labels(lat_labels)
lon_values <- parse_numeric_labels(lon_labels)

# Build results table
results <- data.frame(
  hit_number = seq_len(nrow(hits)),
  row_index = hits[, "row"],
  col_index = hits[, "col"],
  land_mask_value = land_mat[hits],
  ratio_grid_value = ratio_mat[hits],
  latitude = if (length(lat_values) > 0) lat_values[hits[, "row"]] else NA_real_,
  longitude = if (length(lon_values) > 0) lon_values[hits[, "col"]] else NA_real_
)

cat("\nMatched cell(s) for value", target_value, "in land mask:\n")
print(results)

# Optional: write output for downstream use
output_path <- "matched_grid_cells_for_294_9545.csv"
write.csv(results, output_path, row.names = FALSE)
cat("\nSaved match table to:", output_path, "\n")

# Helpful note for interpretation
if (all(is.na(results$latitude)) || all(is.na(results$longitude))) {
  cat(
    "\nNote: Latitude/longitude were not detected from row/column labels.\n",
    "You still have exact matrix coordinates via row_index and col_index.\n",
    sep = ""
  )
}
