unify_metadata <- function(projects) {

  message("Metadata successfully merged.")

  # Construct the destination folder using the global variable norm_method
  #destination_folder <- paste0("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method)
  destination_folder <- paste0("/content/drive/MyDrive/H2Omics_workshop/sequencing_data/", norm_method)

  # List files that match the metadata TSV pattern
  metadata_file <- list.files(destination_folder, pattern = "metadata_formatted\\.tsv$", full.names = TRUE)

  # Read the TSV file (assumes tab-separated values)
  metadata_tsv <- read.delim(metadata_file, col.names = TRUE, show_col_type = FALSE)

  # Display the first 10 rows of the metadata
  print(head(metadata_tsv, 10))

  # Construct the new folder path (using base_path and projects) for input_data
  new_folder <- file.path(base_path, projects, "input_data")

  # Define the new file path for the metadata file (e.g., "ProjectName_metadata_formatted.tsv")
  new_file_path <- file.path(new_folder, paste0(projects, "_metadata_formatted.tsv"))

  # Copy the metadata file to the new folder, overwriting if necessary
  file.copy(from = metadata_file, to = new_file_path, overwrite = TRUE)
}
