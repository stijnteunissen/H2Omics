#' @export
decontam =  function(physeq = resolved_tree_physeq,
                     decon_method = c("frequency", "prevalence", "both"),
                     blank = TRUE) {

  log_message(paste("Step 5: Decontam: Removing contamination.", paste(projects, collapse = ", ")), log_file)

  # Construct the destination folder using the global variable norm_method
  #destination_folder <- paste0("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method)
  destination_folder <- paste0("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data/", norm_method)

  # List files that match the phyloseq RDS pattern for resolved trees
  phyloseq_file <- list.files(destination_folder, pattern = "phyloseq_asv_level_decontam\\.rds$", full.names = TRUE)

  # Read the existing phyloseq file
  phyloseq <- readRDS(phyloseq_file)

  # Construct the new folder path for Before_cleaning_rds_files
  new_folder <- file.path(base_path, projects, "output_data/rds_files/Before_cleaning_rds_files")

  # Define the new file path for saving the resolved phyloseq object as an RDS file
  new_file_path <- file.path(new_folder, paste0(projects, "_phyloseq_asv_level_decontam.rds"))

  # Save the resolved phyloseq object as an RDS file
  saveRDS(phyloseq, file = new_file_path)

  message("Decontam was successfully executed.")

  log_message("Decontam successfully executed.", log_file)
}
