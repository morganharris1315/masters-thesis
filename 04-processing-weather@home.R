# -------------------------------------------------------------------------
# 04-processing-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# Loading in .nc files and understanding file structure
# Calculating RX1day, 33rd Percentile and exceedance days 
# Getting the proportion of greater than or equal to 4 exceedance days
# Creating probability ratio (Future projection/ Current Day) 
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

# 4 dimiensions: longitude0 (44), latitude0 (44), z0 (1), time1 (360)
# 8 vaiables, with rain in `item5216_daily_mean'

var.get.nc(nc,"longitude0")
var.get.nc(nc,"latitude0")

# Have rotated-model grid coordinates (not plain lat/lon)
# longitude0 runs ~194.52 to 213.44, latitude0 runs ~-4.84 to -23.76.

var.inq.nc(nc, "item5216_daily_mean")
precip <- var.get.nc(nc, "item5216_daily_mean")
dim(precip)

# `item5216_daily_mean` is precipitation with dimensions
# Dimensions are (longitude0, latitude0, z0, time1). But as z0 is length 1, R drops it when read and dimesions are...
# 44 x 44 x 360

#Random single grid cell
single_gird_cell <- precip[10, 20, ]

plot(single_gird_cell, type = "l")
max(single_gird_cell, na.rm = TRUE)

#Getting RX1day from precip varaible 
rx1day <- apply(precip, c(1,2), max, na.rm = TRUE)
dim(rx1day)

time <- var.get.nc(nc, "time1")
head(time)
tail(time)
# Time uses a 360_day calendar and units "days since 1840-12-01 00:00:00".
# Values ending in .5 are expected from daily means centred on each day.

close.nc(nc)

# Calculating RX1day -------------------------------------------------------
# For all years (files), all grid cells and both current day and future projections

# All current day files
current_day_files <- list.files ("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/current_decade",
                                 pattern = "\\.nc$",
                                 full.names = TRUE)

# Function to caculate RX1day
compute_rx1day_NetCDF <- function(file) {
  nc <- open.nc(file)
  pr <- var.get.nc(nc, "item5216_daily_mean")
    rx <- apply(pr, c(1,2), max, na.rm = TRUE)
  close.nc(nc)
  return(rx)
}

#Applying function to current day files
current_rx_list <- lapply(current_day_files, compute_rx1day_NetCDF)

length(current_rx_list)
# the length is 3226 e.g there is an RX1day value for every one of the 3226 years. 

current_rx_array <- simplify2array(current_rx_list)
dim(current_rx_array)
# 44 X 44 X 3226 

future_day_files <- list.files("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/model_data/3k_warmer",
                    pattern = "\\.nc$",
                    full.names = TRUE)

future_rx_list <- lapply(future_day_files, compute_rx1day_NetCDF)
length(future_rx_list)
# the length is 2535 e.g there is an RX1day value for every one of the 2535 years. 

future_rx_array <- simplify2array(future_rx_list)
dim(future_rx_array)
# 44 X 44 X 2535

