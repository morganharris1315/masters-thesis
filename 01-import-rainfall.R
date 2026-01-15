# -------------------------------------------------------------------------
# 01-import-rainfall.R 
# -------------------------------------------------------------------------

# helper: extract station number (agent_no) from filename -----------------
extract_agent_no <- function(file_path) {
  file_name <- path_file(file_path)
  agent_no <- str_extract(file_name, "^[0-9]+")
  return(agent_no)
}


# helper: clean region names for file output ------------------------------
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
      col_types = cols(.default = "c"),  # read all as character
      .name_repair = "minimal"
    ) %>%
      clean_names() %>%
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
  
  # Convert key numeric columns
  region_joined <- region_joined %>%
    mutate(
      rainfall_mm = as.numeric(rainfall_mm),
      period_hrs   = as.numeric(period_hrs),
      deficit_mm   = as.numeric(deficit_mm),
      runoff_mm    = as.numeric(runoff_mm),
      latitude     = as.numeric(latitude),
      longitude    = as.numeric(longitude),
      start_year   = as.numeric(start_year),
      end_year     = as.numeric(end_year),
      record_length = as.numeric(record_length)
    )
  
  return(region_joined)
}


# process all regions & save outputs --------------------------------------
process_all_regions <- function(base_raw_dir) {
  
  # Load metadata and clean
  sites_path <- glue("{base_raw_dir}/Identified_Sites.csv")
  message("Loading metadata: ", sites_path)
  
  sites_df <- read_csv(sites_path, show_col_types = FALSE) %>%
    clean_names() %>%
    # Remove any unnamed columns automatically
    select(-starts_with("x")) %>%
    mutate(agent_no = as.character(agent_no))  # fix type mismatch
  
  # Define regions
  region_list <- tribble(
    ~region_name,       ~region_dir,
    "Coromandel",       glue("{base_raw_dir}/coromandel_raw"),
    "Far North",        glue("{base_raw_dir}/far_north_raw"),
    "Top of South",     glue("{base_raw_dir}/top_of_south_raw"),
    "Waikato",          glue("{base_raw_dir}/waikato_raw")
  )
  
  # Output folder
  processed_dir <- glue("{base_raw_dir}/processed")
  dir_create(processed_dir)
  
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
    out_path <- glue("{processed_dir}/rain_{r}_clean.csv")
    
    if (!is.null(out_df)) {
      write_csv(out_df, out_path)
      message("Saved cleaned region file: ", out_path)
    }
  }
  
  # Combine all regions
  all_regions <- bind_rows(region_outputs)
  
  # Save combined dataset
  final_path <- glue("{processed_dir}/rain_all_regions_clean.csv")
  write_csv(all_regions, final_path)
  
  message("\n✔ Complete. All rainfall datasets processed successfully.")
  
  return(all_regions)
}


# Run with: ---------------------------------------------------------------
process_all_regions(base_raw_dir)


