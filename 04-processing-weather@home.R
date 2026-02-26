# -------------------------------------------------------------------------
# 04-processing-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Loading in .nc files
# Understanding file structure
# Caculating RX1day
# Caculating 33rd Percentile 
# Caculating Exceedance days 
# Proportion of greater than or equal to 4 exceedance days
# Probabilty ratio (Future projection/ Current Day) 
# -------------------------------------------------------------------------


# Loading in RNetCDF ------------------------------------------------------
#install.packages("RNetCDF")
library(RNetCDF)

# Exploring the structure of the data -------------------------------------
#Looking at a single file in the list of files in the folder current_day
current_day_files <- list.files ("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/current_decade",
                      pattern = "\\.nc$",
                      full.names = TRUE)

nc <- open.nc(current_day_files[1])

file.inq.nc(nc)

file.inq.nc(nc)['format']

print.nc(nc)

var.get.nc(nc,"longitude0")
var.get.nc(nc,"latitude0")

var.inq.nc(nc, "item5216_daily_mean")
pr <- var.get.nc(nc, "item5216_daily_mean")

dim(pr)

single_gird_cell <- pr[10, 20, ]

plot(single_gird_cell, type = "l")
max(single_gird_cell, na.rm = TRUE)

rx1day <- apply(pr, c(1,2), max, na.rm = TRUE)

dim(rx1day)

time <- var.get.nc(nc, "time1")
head(time)
tail(time)

close.nc(nc)


# Calculating RX1day -------------------------------------------------------
# For all years (files), all grid cells and both current day and future projections

current_day_files <- list.files ("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/current_decade",
                                 pattern = "\\.nc$",
                                 full.names = TRUE)

compute_rx1day_NetCDF <- function(file) {
  nc <- open.nc(file)
  
  pr <- var.get.nc(nc, "item5216_daily_mean")
  
  # pr is already [44,44,360]
  rx <- apply(pr, c(1,2), max, na.rm = TRUE)
  
  close.nc(nc)
  
  return(rx)
}

current_rx_list <- lapply(current_day_files, compute_rx1day_NetCDF)

length(current_rx_list)

current_rx_array <- simplify2array(current_rx_list)
dim(current_rx_array)

future_day_files <- list.files(
  "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/3k_warmer",
  pattern = "\\.nc$",
  full.names = TRUE
)

future_rx_list <- lapply(future_day_files, compute_rx1day_NetCDF)
future_rx_array <- simplify2array(future_rx_list)

dim(future_rx_array)
