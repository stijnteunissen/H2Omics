create_folders = function(projects) {

  for (project in projects) {

    project_folder = paste0(base_path, project)

    # Create the required directories if they don't exist
    dir.create(paste0(project_folder, "/input_data"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(project_folder, "/output_data"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(project_folder, "/output_data/csv_files"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(project_folder, "/output_data/rds_files"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(project_folder, "/messages"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(project_folder, "/figures"), showWarnings = FALSE, recursive = TRUE)

    # Define source (git_folder) and destination (qiime2_output) for the initial copy
    #git_folder <- paste0("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method)
    git_folder <- paste0("/content/drive/MyDrive/H2Omics_workshop/sequencing_data/", norm_method)
    source_folder <- paste0(project_folder, "/qiime2_output")
    destination_folder <- paste0(project_folder, "/input_data")

    # List files in the git folder
    git_files = list.files(git_folder, full.names = TRUE)

    # Define the patterns of the files that need to be copied
    required_patterns = c("table.*\\.qza$", "rooted-tree.*\\.qza$", "classifier.*\\.qza$", "metadata\\.tsv$")

    # Filter files matching the required patterns
    files_to_copy = git_files[sapply(required_patterns, function(pattern) any(grepl(pattern, git_files)))]

    # Copy relevant files from git_folder to qiime2_output
    file.copy(files_to_copy, source_folder, overwrite = TRUE)

    # List all files in source folder again (after copying)
    files = list.files(source_folder, full.names = TRUE)

    # Check for required files - stop if missing
    for (file_pattern in required_patterns) {
      if (!any(grepl(file_pattern, files))) {
        error_message = paste("Error:", file_pattern, "does not exist in", source_folder, "for project:", project, "\n")
        log_message(error_message, log_file)
        stop(error_message)
      }
    }

    # Check for optional files - only warning if missing
    optional_files = c("qPCR.*\\.csv$", "fcm.*\\.csv$", "metadata_extra\\.tsv$", "prediction*\\.RDS$")
    for (file_pattern in optional_files) {
      if (!any(grepl(file_pattern, files))) {
        warning_message = paste("Warning:", file_pattern, "does not exist in", source_folder, "for project:", project, "\n")
        log_message(warning_message, log_file)
      }
    }

    # Copy relevant files from qiime2_output to input_data
    files_for_phyloseq_object = list.files(
      source_folder,
      pattern = "table.*\\.qza$|rooted-tree.*\\.qza$|classifier.*\\.qza$|metadata\\.tsv$|metadata_extra\\.tsv$|dna-sequences.*\\.csv$|fcm.*\\.csv$|qPCR.*\\.csv$|prediction*\\.RDS$",
      full.names = TRUE
    )

    file.copy(files_for_phyloseq_object, destination_folder, overwrite = TRUE)

  }
}
