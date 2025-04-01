rarefying <- function(physeq = physeq,
                      norm_method = NULL,
                      iteration = 100) {

  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder using the global variable norm_method
  destination_folder <- file.path("/content/drive/MyDrive/H2Omics_workshop/sequencing_data", norm_method, cc_val)

  # Define the source folder (After_cleaning_rds_files)
  source_folder <- file.path(destination_folder, "After_cleaning_rds_files")

  # Define the destination folder where the entire directory should be copied
  target_folder <- file.path(base_path, projects, "output_data", "After_cleaning_rds_files")

  # Check if the source folder exists
  if (dir.exists(source_folder)) {
  # Copy the entire directory and its contents to the target location
  file.copy(from = source_folder, to = target_folder, recursive = TRUE, overwrite = TRUE)
  }

  message("Data has been rarefied.")
}
