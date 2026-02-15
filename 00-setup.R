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
model_data_dir       <- glue("{base_raw_dir}/model_data")

theme_thesis <- theme_minimal(base_size = 9) +
  theme(
    plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 9),
    axis.text    = element_text(size = 7),
    strip.text   = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    panel.spacing = unit(0.5, "lines")
  )

fig_width_full  <- 6.7
fig_height_short <- 3.8
fig_height_med   <- 4.8
fig_height_tall  <- 6.2
fig_height_full <- 7

fig_width_standard  <- 16.89 / 2.54
fig_height_standard <- 19 / 2.54

fig_width_hoz_full  <- 11
fig_height_hoz_full <- 7.2

fig_width_hoz_half  <- 8.5
fig_height_hoz_half <- 6

