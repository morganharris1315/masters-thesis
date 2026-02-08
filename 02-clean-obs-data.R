# -------------------------------------------------------------------------
# 02-clean-obs-data.R
# -------------------------------------------------------------------------

processed_dir_obs <- glue("{base_raw_dir}/obs_data/cleaned_datasets")

raw_obs_dirs <- list(
  coromandel   = coromandel_raw_dir,
  far_north    = far_north_raw_dir,
  top_of_south = top_of_south_raw_dir,
  waikato      = waikato_raw_dir
)

regions_obs <- c("coromandel", "far_north", "top_of_south", "waikato")

# Columns expected to be numeric
numeric_cols <- c("rainfall_mm", "period_hrs", "deficit_mm", "runoff_mm",
                  "latitude", "longitude", "start_year", "end_year", "record_length")


# Load Identified Sites metadata ------------------------------------------
sites_path <- glue("{base_raw_dir}/Identified_Sites.csv")

sites_df <- read_csv(sites_path, show_col_types = FALSE) %>%
  clean_names() %>%
  select(-starts_with("x")) %>%  # remove unnamed columns
  mutate(agent_no = as.character(agent_no))


# function: extract agent_no from filename --------------------------------
extract_agent_no <- function(file_path) {
  fname <- path_file(file_path)
  str_extract(fname, "^[0-9]+")
}


# Function: clean a single region -----------------------------------------
clean_region_timeseries_obs <- function(region_name) {
  
  message(glue("Processing region: {region_name}"))
  
  csv_files <- dir_ls(raw_obs_dirs[[region_name]], glob = "*.csv")
  if(length(csv_files) == 0){
    message(glue("No CSV files found for region {region_name}"))
    return(NULL)
  }
  
  # Read all files as character to avoid bind_rows type conflicts
  df_list <- map(csv_files, function(f) {
    agent_no <- extract_agent_no(f)
    
    read_csv(f,
             col_types = cols(.default = "c"),
             show_col_types = FALSE) %>%
      clean_names() %>%
      mutate(agent_no = agent_no)
  })
  
  df <- bind_rows(df_list)
  
  # Join with sites_df to get station names
  df <- df %>%
    left_join(sites_df, by = "agent_no")
  
  # Ensure 'station' column exists
  if(!"station" %in% colnames(df)){
    df$station <- paste0("unknown_station")
  }
  
  # Standardize timestamp column
  timestamp_col <- intersect(c("observation_time_utc", "observation_time", "datetime", "date_time"),
                             colnames(df))
  
  if(length(timestamp_col) == 0){
    message(glue("⚠ Region {region_name} has no timestamp column. Creating observation_date = NA"))
    df$observation_date <- NA
  } else {
    ts_col <- timestamp_col[1]
    df <- df %>%
      mutate(
        observation_time_utc = as.character(.data[[ts_col]]),
        observation_date = as.Date(substr(.data[[ts_col]], 1, 10))
      ) %>%
      filter(!is.na(observation_date))
  }
  
  # Convert numeric columns
  for(col in numeric_cols){
    if(col %in% colnames(df)) df[[col]] <- as.numeric(df[[col]])
  }
  
  # Convert known mixed-type metadata columns to character
  mixed_cols <- c("case_study", "data_source_water_balance")
  for(col in mixed_cols){
    if(col %in% colnames(df)) df[[col]] <- as.character(df[[col]])
  }
  
  # Fill missing dates per station
  df_filled <- df %>%
    group_by(station) %>%
    complete(observation_date = seq.Date(min(observation_date), max(observation_date), by = "day")) %>%
    ungroup() %>%
    mutate(region = region_name)
  
  # Fill numeric columns for added rows
  for(col in numeric_cols){
    if(col %in% colnames(df_filled)) df_filled[[col]][is.na(df_filled[[col]])] <- NA
  }
  
  # Missingness summary
  message(glue("Missingness summary for {region_name}:"))
  missing_summary <- df_filled %>%
    group_by(station) %>%
    summarise(
      first_date = min(observation_date, na.rm = TRUE),
      last_date = max(observation_date, na.rm = TRUE),
      total_days = n(),
      missing_rainfall = sum(is.na(rainfall_mm)),
      .groups = "drop"
    )
  print(missing_summary)
  
  # Save cleaned CSV
  output_file <- glue("{processed_dir_obs}/rain_{region_name}_cleaned.csv")
  write_csv(df_filled, output_file)
  message(glue("Saved cleaned CSV for region {region_name}: {output_file}"))
  
  return(df_filled)
}


# process all regions -----------------------------------------------------
all_regions_obs <- map(regions_obs, clean_region_timeseries_obs)
names(all_regions_obs) <- regions_obs

# Combine all regions
combined_df_obs <- bind_rows(all_regions_obs, .id = "region")

# Add hydrological year (July → June)
combined_df_obs <- combined_df_obs %>%
  mutate(hydro_year = if_else(month(observation_date) >= 7,
                              year(observation_date) + 1,
                              year(observation_date)))

# Missingness report
missingness_report <- combined_df_obs %>%
  group_by(region, station) %>%
  summarise(
    first_date = min(observation_date, na.rm = TRUE),
    last_date = max(observation_date, na.rm = TRUE),
    total_days = n(),
    missing_rainfall = sum(is.na(rainfall_mm)),
    pct_missing = round(100 * missing_rainfall / total_days, 2),
    .groups = "drop"
  ) %>%
  arrange(region, station)

print(missingness_report, n = 50)

# Save final combined dataset
write_csv(combined_df_obs, glue("{processed_dir_obs}/rain_all_regions_cleaned.csv"))


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

# identifying the best stations
least_missing_stations <- missingness_report %>%
  arrange(pct_missing) %>%
  print(least_missing_stations, n = 40, width = Inf)

longest_record <- missingness_report %>%
  arrange(desc(total_days)) %>%
  print(longest_record, n = 40, width = Inf)

# check stations in each region
unique(combined_df_obs$station[combined_df_obs$region == "coromandel"])
unique(combined_df_obs$station[combined_df_obs$region == "far_north"])
unique(combined_df_obs$station[combined_df_obs$region == "top_of_south"])
unique(combined_df_obs$station[combined_df_obs$region == "waikato"])
