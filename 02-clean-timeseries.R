# -------------------------------------------------------------------------
# 02-clean-timeseries.R
# -------------------------------------------------------------------------

# processed directory -----------------------------------------------------
processed_dir <- glue("{base_raw_dir}/processed")

regions <- c("coromandel", "far_north", "top_of_south", "waikato")


# Columns that should be numeric ------------------------------------------
numeric_cols <- c("rainfall_mm", "period_hrs", "deficit_mm", "runoff_mm",
                  "latitude", "longitude", "start_year", "end_year", "record_length")


# Function to clean a single region ---------------------------------------
clean_region_timeseries <- function(region_name) {
  
  message(glue("Processing region: {region_name}"))
  
  # Path to CSV
  csv_file <- glue("{processed_dir}/rain_{region_name}_clean.csv")
  
  if (!file_exists(csv_file)) {
    message(glue("No CSV file found for region: {region_name}"))
    return(NULL)
  }
  
  df <- read_csv(csv_file, show_col_types = FALSE)
  
  # Convert observation_time_utc to Date
  df <- df %>%
    mutate(
      observation_time_utc = as.character(observation_time_utc),
      observation_date = as.Date(substr(observation_time_utc, 1, 10))
    ) %>%
    filter(!is.na(observation_date))  # Remove invalid rows if any
  
  if (nrow(df) == 0) {
    message(glue("Warning: Region {region_name} has no valid observation dates"))
    return(NULL)
  }
  
  # numeric columns
  for (col in numeric_cols) {
    if (col %in% colnames(df)) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }
  
  # Fill missing dates per station
  df_filled <- df %>%
    group_by(station) %>%
    complete(observation_date = seq.Date(min(observation_date), 
                                         max(observation_date), 
                                         by = "day")) %>%
    ungroup() %>%
    mutate(
      region = region_name
    )
  
  # Fill numeric columns with NA for added rows
  for (col in numeric_cols) {
    if (col %in% colnames(df_filled)) {
      df_filled[[col]][is.na(df_filled[[col]])] <- NA
    }
  }
  
  # Missingness summary
  message(glue("Missingness summary for {region_name}:"))
  missing_summary <- df_filled %>%
    group_by(station) %>%
    summarise(
      first_date = min(observation_date, na.rm = TRUE),
      last_date  = max(observation_date, na.rm = TRUE),
      total_days = n(),
      missing_rainfall = sum(is.na(rainfall_mm)),
      .groups = "drop"
    )
  print(missing_summary)
  
  # Save cleaned file
  output_file <- glue("{processed_dir}/rain_{region_name}_cleaned.csv")
  write_csv(df_filled, output_file)
  message(glue("Saved cleaned CSV for region {region_name}: {output_file}"))
  
  return(df_filled)
}

# process all regions -----------------------------------------------------
all_regions <- map(regions, clean_region_timeseries)
names(all_regions) <- regions

# combine all regions 
combined_df <- bind_rows(all_regions, .id = "region")

# creating hydrological year ----------------------------------------------
combined_df <- combined_df %>%
  mutate(
    hydro_year = if_else(
      month(observation_date) >= 7,   # July → December
      year(observation_date) + 1,     # Belongs to next hydrological year
      year(observation_date)          # Jan → June belongs to current year
    )
  )

# combined missingness report ---------------------------------------------
missingness_report <- combined_df %>%
  group_by(region, station) %>%
  summarise(
    first_date = min(observation_date, na.rm = TRUE),
    last_date  = max(observation_date, na.rm = TRUE),
    total_days = n(),
    missing_rainfall = sum(is.na(rainfall_mm)),
    pct_missing = round(100 * missing_rainfall / total_days, 2),
    .groups = "drop"
  ) %>%
  arrange(region, station)

# Print top rows to console
print(missingness_report, n = 50)

# top stations by percentage of missing rainfall
top_missing_stations <- missingness_report %>%
  arrange(desc(pct_missing)) %>%
  slice_head(n = 20)  #adjust n
print(top_missing_stations)

# stations with more than 20% missing rainfall
over20_missing_stations <- missingness_report %>%
  filter(pct_missing > 20) %>%
  arrange(desc(pct_missing))
print(over20_missing_stations)
#there are seven stations with over 20% missing (Parengarenga, Mangarata, Matamata,Calders Fm, Kerikeri Aero 2, Kohatu,Highfield 2 and Stony Bay)

#identifying the best stations
least_missing_stations <- missingness_report %>%
  arrange(pct_missing) %>%
print(least_missing_stations, n = 40, width = Inf)

longest_record <- missingness_report %>%
  arrange(desc(total_days)) %>%
  print(longest_record, n = 40, width = Inf)


unique(combined_df$station[combined_df$region == "coromandel"])

