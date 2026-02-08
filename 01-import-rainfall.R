# -------------------------------------------------------------------------
# 01-import-obs-data.R
# -------------------------------------------------------------------------

# helper: extract agent number from filename ------------------------------
extract_agent_no <- function(file_path) {
  file_name <- path_file(file_path)
  agent_no <- str_extract(file_name, "^[0-9]+")
  return(agent_no)
}

# helper: clean region names for output ----------------------------------
clean_region_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

# function: read all rainfall CSV files in a region directory -------------
read_region_rainfall <- function(region_dir, region_name, sites_df) {
  
  message("\n------------------------------------------------------------")
  message("Reading region: ", region_name)
  
  files <- dir_ls(region_dir, regexp = "\\.csv$")
  message("Found ", length(files), " rainfall files.")
  
  if (length(files) == 0) {
    warning("No files found for region: ", region_name)
    return(NULL)
  }
  
  station_dfs <- map(files, function(file) {
    
    agent_no <- extract_agent_no(file)
    message("  • Reading station ", agent_no)
    
    df <- vroom(
      file,
      show_col_types = FALSE,
      col_types = cols(.default = "c"),
      .name_repair = "minimal"
    ) %>%
      janitor::clean_names() %>%
      mutate(
        agent_no = agent_no,
        region = region_name
      )
    
    df
  })
  
  # Combine all stations for this region
  region_df <- bind_rows(station_dfs)
  message("Combined ", length(station_dfs), " station files for ", region_name)
  
  # Join with metadata
  region_joined <- region_df %>%
    left_join(sites_df, by = "agent_no")
  
  message("Joined rainfall data with Identified_Sites metadata.")
  
  # Convert key numeric columns safely
  numeric_cols <- c("rainfall_mm", "period_hrs", "deficit_mm", "runoff_mm",
                    "latitude", "longitude", "start_year", "end_year", "record_length")
  
  for (col in numeric_cols) {
    if (col %in% colnames(region_joined)) {
      region_joined[[col]] <- as.numeric(region_joined[[col]])
    }
  }
  
  return(region_joined)
}

# Process all regions & save outputs --------------------------------------
process_all_regions <- function(base_raw_dir) {
  
  # Load metadata
  sites_path <- glue("{base_raw_dir}/Identified_Sites.csv")
  message("Loading metadata: ", sites_path)
  
  sites_df <- read_csv(sites_path, show_col_types = FALSE) %>%
    clean_names() %>%
    select(-starts_with("x")) %>%
    mutate(agent_no = as.character(agent_no))
  
  # Define regions
  region_list <- tribble(
    ~region_name,       ~region_dir,
    "Coromandel",       glue("{base_raw_dir}/coromandel_raw"),
    "Far North",        glue("{base_raw_dir}/far_north_raw"),
    "Top of South",     glue("{base_raw_dir}/top_of_south_raw"),
    "Waikato",          glue("{base_raw_dir}/waikato_raw")
  )
  
  # Output folder
  processed_dir_obs <- glue("{base_raw_dir}/obs_data/cleaned_datasets")
  
  # Process each region
  region_outputs <- map2(
    region_list$region_dir,
    region_list$region_name,
    ~ read_region_rainfall(.x, .y, sites_df)
  )
  
  names(region_outputs) <- map_chr(region_list$region_name, clean_region_name)
  
  # Save region-level files
  for (r in names(region_outputs)) {
    out_df <- region_outputs[[r]]
    out_path <- glue("{processed_dir_obs}/rain_{r}_clean.csv")
    
    if (!is.null(out_df)) {
      write_csv(out_df, out_path)
      message("Saved cleaned region file: ", out_path)
    }
  }
  
  # Combine all regions
  all_regions <- bind_rows(region_outputs)
  
  # Save combined dataset
  final_path <- glue("{processed_dir_obs}/rain_all_regions_clean.csv")
  write_csv(all_regions, final_path)
  
  message("\n Complete. All rainfall datasets processed successfully.")
  
  return(all_regions)
}

# Run ---------------------------------------------------------------------
all_regions <- process_all_regions(base_raw_dir)




