#' @export
rarefying <- function(physeq = physeq,
                      norm_method = NULL,
                      iteration = 100) {

  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder using the global variable norm_method
  destination_folder <- file.path("/content/Workshop_H2Omics_test/H2Omics_workshopp/sequencing_data", norm_method, cc_val)

  # Define the source folder (After_cleaning_rds_files)
  source_folder <- file.path(destination_folder, "After_cleaning_rds_files")

  # Define the destination folder where the entire directory should be copied
  target_folder <- file.path(base_path, projects, "output_data", "rds_files", "After_cleaning_rds_files")

  # Check if the source folder exists
  if (dir.exists(source_folder)) {
  # Copy the entire directory and its contents to the target location
  file.copy(from = source_folder, to = target_folder, recursive = TRUE, overwrite = TRUE)
  }

  # Define the source folder (After_cleaning_rds_files)
  source_folder_csv <- file.path(destination_folder, "csv_files")

  # Define the destination folder where the entire directory should be copied
  target_folder_csv <- file.path(base_path, projects, "output_data", "csv_files")

  # Check if the source folder exists
  if (dir.exists(source_folder_csv)) {
    # Copy the entire directory and its contents to the target location
    file.copy(from = source_folder_csv, to = target_folder_csv, recursive = TRUE, overwrite = TRUE)
  }

  message("Data has been rarefied.")
}
