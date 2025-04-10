#' @export
creating_physeq_object = function(projects) {

  message("Phyloseq successfully created.")

  # Construct the destination folder using the global variable norm_method
  #destination_folder <- paste0("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method)
  destination_folder <- paste0("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data/", norm_method)

  # List files that match the phyloseq RDS pattern
  phyloseq_file <- list.files(destination_folder, pattern = "phyloseq_uncleaned\\.rds$", full.names = TRUE)

  # Read the RDS file containing the phyloseq object
  phyloseq <- readRDS(phyloseq_file)

  # Construct the new folder path (using base_path and projects) for Before_cleaning_rds_files
  new_folder <- file.path(base_path, projects, "output_data/rds_files/Before_cleaning_rds_files")

  # Define the new file path for saving the phyloseq object as an RDS file
  new_file_path <- file.path(new_folder, paste0(projects, "_phyloseq_uncleaned.rds"))

  # Save the phyloseq object as an RDS file
  saveRDS(phyloseq, file = new_file_path)

  print(phyloseq)
}
