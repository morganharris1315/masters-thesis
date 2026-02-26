# -------------------------------------------------------------------------
# 04-processing-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# 
# -------------------------------------------------------------------------

if (!requireNamespace("RNetCDF", quietly = TRUE)) {
  install.packages("RNetCDF")
}
library(RNetCDF)

current_decade_CDF <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/obs_data/current_decade.nc"

if (!file.exists(current_decade_CDF)) {
  stop(
    sprintf(
      "NetCDF file not found: %s\nCheck the path and filename (including the .nc extension).",
      current_decade_CDF
    )
  )
}

if (file.access(current_decade_CDF, mode = 4) != 0) {
  stop(
    sprintf(
      "No read permission for file: %s\nMove it out of a restricted folder (for example OneDrive) or update file permissions.",
      current_decade_CDF
    )
  )
}

current_decade <- open.nc(current_decade_CDF)

file.inq.nc(current_decade)

file.inq.nc(current_decade)['format']

print.nc(current_decade)
