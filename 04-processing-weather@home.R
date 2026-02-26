# -------------------------------------------------------------------------
# 04-processing-weather@home.R
# -------------------------------------------------------------------------
# Feb 2026
# 
# -------------------------------------------------------------------------

install.packages("RNetCDF")
library(RNetCDF)

current_decade_CDF <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Compound Events/obs_data/current_decade"
current_decade <- open.nc(current_decade_CDF)

file.inq.nc(current_decade)

file.inq.nc(current_decade)['format']

print.nc(current_decade)
