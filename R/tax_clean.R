#' @export
tax_clean = function(physeq = physeq, tax_filter = TRUE) {

  # Construct the destination folder using the global variable norm_method
  #destination_folder <- paste0("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method)
  destination_folder <- paste0("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data/", norm_method)

  # List files that match the phyloseq RDS pattern
  uncleaned_phyloseq_file <- list.files(destination_folder, pattern = "phyloseq_uncleaned\\.rds$", full.names = TRUE)

  # List files that match the phyloseq cleaned RDS pattern
  cleaned_phyloseq_file <- list.files(destination_folder, pattern = "phyloseq_cleaned\\.rds$", full.names = TRUE)

  # Read the RDS file containing the phyloseq object
  uncleaned_phyloseq <- readRDS(uncleaned_phyloseq_file)
  cleaned_phyloseq <- readRDS(cleaned_phyloseq_file)

  print(head(tax_table(uncleaned_phyloseq), 10))
  print(head(tax_table(cleaned_phyloseq), 10))


  # # Convert OTU tables to matrices for subsetting
  # otu_mat_uncleaned <- as(tax_tabel(uncleaned_phyloseq), "matrix")
  # otu_mat_cleaned <- as(tax_tabel(cleaned_phyloseq), "matrix")
  #
  # # Determine subset dimensions (first 10 OTUs and first 5 samples)
  # n_taxa_uncleaned <- min(10, nrow(otu_mat_uncleaned))
  # n_samples_uncleaned <- min(4, ncol(otu_mat_uncleaned))
  #
  # n_taxa_cleaned <- min(10, nrow(otu_mat_cleaned))
  # n_samples_cleaned <- min(4, ncol(otu_mat_cleaned))
  #
  # # Print the subsets
  # cat("First 10 OTUs and first 5 samples of uncleaned phyloseq object:\n")
  # print(otu_mat_uncleaned[1:n_taxa_uncleaned, 1:n_samples_uncleaned])
  #
  # cat("\nFirst 10 OTUs and first 5 samples of cleaned phyloseq object:\n")
  # print(otu_mat_cleaned[1:n_taxa_cleaned, 1:n_samples_cleaned])

  # Construct the new folder path (using base_path and projects) for Before_cleaning_rds_files
  new_folder <- file.path(base_path, projects, "output_data/rds_files/Before_cleaning_rds_files")

  # Define the new file path for saving the cleaned phyloseq object as an RDS file
  new_file_path <- file.path(new_folder, paste0(projects, "_phyloseq_cleaned.rds"))

  # Save the cleaned phyloseq object as an RDS file
  saveRDS(cleaned_phyloseq, file = new_file_path)

  # Print a message indicating the process is complete
  message("Phyloseq object has been cleaned and saved.")
}
