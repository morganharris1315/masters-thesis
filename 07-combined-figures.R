# -------------------------------------------------------------------------
# 07-combined-figures.R
# -------------------------------------------------------------------------
# Mar 2026
# Find the weather@home grid indices (lon_index/lat_index) for a specific
# land-sea mask value. Uses the same mask transpose logic as
# 05-mapping-weather@home.R.
# -------------------------------------------------------------------------

# Setting input paths and target value ------------------------------------
model_data_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data"
land_mask_file <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/Land-Sea Mask for Weather@home Data.csv"
ratio_grid_file <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/weather@home_exceedance_ge4_ge5_top10_joint_probability_ratio_grid.csv"
nc_file <- file.path(
  model_data_dir,
  "current_decade",
  "item5216_daily_mean_a000_2006-07_2007-06-NZtrim-mm.nc"
)

mask_transpose <- TRUE
target_value <- 294.9544983

# Helper functions ---------------------------------------------------------
find_matching_cells_exact <- function(grid_matrix, value) {
  which(grid_matrix == value, arr.ind = TRUE)
}

load_first_nc_slice <- function(file_path, variable_name = "item5216_daily_mean") {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required to read the NetCDF file.")
  }
  
  nc <- ncdf4::nc_open(file_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  var_info <- nc$var[[variable_name]]
  if (is.null(var_info)) {
    stop("Variable '", variable_name, "' was not found in NetCDF file: ", file_path)
  }
  
  var_ndims <- length(var_info$dim)
  if (var_ndims < 2) {
    stop("Expected at least 2 dimensions in NetCDF variable: ", variable_name)
  }
  
  start <- rep(1, var_ndims)
  count <- rep(1, var_ndims)
  count[1:2] <- -1
  
  nc_slice <- ncdf4::ncvar_get(nc, variable_name, start = start, count = count)
  as.matrix(nc_slice)
}

load_lse_mask_matrix <- function(mask_file, nlon, nlat, transpose_mask = FALSE) {
  raw <- utils::read.csv(
    file = mask_file,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NaN", "NA", "")
  )
  
  numeric_raw <- suppressWarnings(as.data.frame(lapply(raw, as.numeric)))
  
  as_m <- function(x) as.matrix(x)
  base_candidates <- list(
    as_m(numeric_raw),
    as_m(numeric_raw[-1, , drop = FALSE]),
    as_m(numeric_raw[, -1, drop = FALSE]),
    as_m(numeric_raw[-1, -1, drop = FALSE])
  )
  
  dims_ok <- vapply(base_candidates, function(m) {
    is.matrix(m) && nrow(m) == nlon && ncol(m) == nlat
  }, logical(1))
  
  if (!any(dims_ok)) {
    stop(sprintf("Could not coerce LSE mask to %d x %d matrix.", nlon, nlat))
  }
  
  mask <- base_candidates[[which(dims_ok)[1]]]
  if (isTRUE(transpose_mask)) {
    mask <- t(mask)
  }
  
  if (nrow(mask) != nlon || ncol(mask) != nlat) {
    stop("Mask dimensions do not match rainfall grid after transpose setting.")
  }
  
  message(sprintf(
    "Loaded LSE mask as %d x %d matrix (transpose=%s).",
    nrow(mask),
    ncol(mask),
    ifelse(isTRUE(transpose_mask), "TRUE", "FALSE")
  ))
  mask
}

# Loading grid files -------------------------------------------------------
nc_grid_matrix <- load_first_nc_slice(nc_file)
nlon <- nrow(nc_grid_matrix)
nlat <- ncol(nc_grid_matrix)
land_mask_matrix <- load_lse_mask_matrix(
  mask_file = land_mask_file,
  nlon = nlon,
  nlat = nlat,
  transpose_mask = mask_transpose
)
ratio_grid <- utils::read.csv(ratio_grid_file, stringsAsFactors = FALSE)

# Checking dimensions ------------------------------------------------------
if (!all(dim(land_mask_matrix) == dim(nc_grid_matrix))) {
  stop(
    "Grid dimensions do not match.\n",
    "Land mask dimensions: ", paste(dim(land_mask_matrix), collapse = " x "), "\n",
    "NetCDF grid dimensions: ", paste(dim(nc_grid_matrix), collapse = " x "), "\n",
    "Check whether one file needs transposing or different row/column handling."
  )
}

if (!all(c("lon_index", "lat_index") %in% names(ratio_grid))) {
  stop("Ratio grid must contain lon_index and lat_index columns.")
}

# Finding matching grid cell(s) -------------------------------------------
match_indices <- find_matching_cells_exact(land_mask_matrix, target_value)

if (nrow(match_indices) == 0) {
  stop(
    "No grid cell found in the land mask for exact value ",
    target_value,
    "."
  )
}

# Building output table ----------------------------------------------------
matched_cells <- data.frame(
  hit_number = seq_len(nrow(match_indices)),
  lon_index = match_indices[, "row"],
  lat_index = match_indices[, "col"],
  land_mask_value = land_mask_matrix[match_indices]
)

ratio_lookup <- ratio_grid[, c("lon_index", "lat_index"), drop = FALSE]
matched_cells$ratio_row_matches <- vapply(seq_len(nrow(matched_cells)), function(i) {
  sum(
    ratio_lookup$lon_index == matched_cells$lon_index[i] &
      ratio_lookup$lat_index == matched_cells$lat_index[i],
    na.rm = TRUE
  )
}, integer(1))

cat("Matching coordinates in console:\n")
print(matched_cells)
