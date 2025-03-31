tax_clean = function(physeq = physeq, tax_filter = TRUE) {

  # Check if BiocManager is installed, and install if necessary
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
  }
  suppressMessages(library(BiocManager))

  # Install phyloseq package if not installed
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    suppressMessages(BiocManager::install("phyloseq", quiet = TRUE, update = TRUE, ask = FALSE))
  }
  suppressMessages(library(phyloseq))

  # Construct the destination folder using the global variable norm_method
  destination_folder <- paste("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method, sep = "")

  # List files that match the phyloseq RDS pattern
  uncleaned_phyloseq_file <- list.files(destination_folder, pattern = "phyloseq_uncleaned\\.RDS$", full.names = TRUE)

  # List files that match the phyloseq cleaned RDS pattern
  cleaned_phyloseq_file <- list.files(destination_folder, pattern = "phyloseq_cleaned\\.RDS$", full.names = TRUE)

  # Read the RDS file containing the phyloseq object
  uncleaned_phyloseq <- readRDS(uncleaned_phyloseq_file)
  cleaned_phyloseq <- readRDS(cleaned_phyloseq_file)

  # Print the first 10 rows of the OTU table for both uncleaned and cleaned phyloseq objects
  print(head(otu_table(uncleaned_phyloseq), 10))
  print(head(otu_table(cleaned_phyloseq), 10))

  # Construct the new folder path (using base_path and projects) for Before_cleaning_rds_files
  new_folder <- file.path(base_path, projects, "output_data/rds_files/Before_cleaning_rds_files")

  # Define the new file path for saving the cleaned phyloseq object as an RDS file
  new_file_path <- file.path(new_folder, paste0(projects, "_phyloseq_cleaned.RDS"))

  # Save the cleaned phyloseq object as an RDS file
  saveRDS(cleaned_phyloseq, file = new_file_path)

  # Print a message indicating the process is complete
  message("Phyloseq object has been cleaned and saved.")
}
