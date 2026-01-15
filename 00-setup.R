# -------------------------------------------------------------------------
# 00-setup.R
# -------------------------------------------------------------------------

#rm(list = ls())

library(tidyverse)
library(readxl)
library(janitor)
library(scales)
library(skimr)
library(patchwork)
library(glue)
library(fs)
library(lubridate)
library(vroom)
library(ggplot2)
library(dplyr)

options(scipen = 999, digits = 4)

# directories --------------------------------------------------------------
thesis_dir <- "C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis"
base_raw_dir <- glue("{thesis_dir}/Historic Compound Events")
coromandel_raw_dir   <- glue("{base_raw_dir}/coromandel_raw")
far_north_raw_dir    <- glue("{base_raw_dir}/far_north_raw")
top_of_south_raw_dir <- glue("{base_raw_dir}/top_of_south_raw")
waikato_raw_dir      <- glue("{base_raw_dir}/waikato_raw")

