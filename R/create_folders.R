create_folders = function(projects) {

  project_folder = paste0(base_path, projects)

    # Create the required directories if they don't exist
    if(!dir.exists(paste0(project_folder, "/input_data"))){dir.create(paste0(project_folder, "/input_data"))}
    if(!dir.exists(paste0(project_folder, "/output_data"))){dir.create(paste0(project_folder, "/output_data"))}
    if(!dir.exists(paste0(project_folder, "/output_data/csv_files"))){dir.create(paste0(project_folder, "/output_data/csv_files"))}
    if(!dir.exists(paste0(project_folder, "/output_data/rds_files"))){dir.create(paste0(project_folder, "/output_data/rds_files"))}
    if(!dir.exists(paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files"))){dir.create(paste0(project_folder, "/output_data/rds_files/Before_cleaning_rds_files"))}
    if(!dir.exists(paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files"))){dir.create(paste0(project_folder, "/output_data/rds_files/After_cleaning_rds_files"))}
    if(!dir.exists(paste0(project_folder, "/messages"))){dir.create(paste0(project_folder, "/messages"))}
    if(!dir.exists(paste0(project_folder, "/figures"))){dir.create(paste0(project_folder, "/figures"))}

    source_folder = paste0(project_folder, "/qiime2_output")
    destination_folder = paste0(project_folder, "/input_data")

    files = list.files(source_folder, full.names = TRUE)

    # Check for required files - stop if missing
    required_files = c("table.*\\.qza$", "rooted-tree.*\\.qza$", "classifier.*\\.qza", "metadata\\.tsv$")
    for (file_pattern in required_files) {
      if (!any(grepl(file_pattern, files))) {  # Use 'files' instead of 'file'
        error_message = paste("Error:", file_pattern, "does not exist in", source_folder, "for project:", projects, "\n")
        log_message(error_message, log_file)
        stop(error_message)
      }
    }

    # Check for optional files - only warning if missing
    optional_files = c("qPCR.*\\.csv$", "fcm.*\\.csv$", "metadata_extra\\.tsv$", "prediction*\\.RDS$")
    for (file_pattern in optional_files) {
      if (!any(grepl(file_pattern, files))) {  # Use 'files' instead of 'file'
        warning_message = paste("Warning:", file_pattern, "does not exist in", source_folder, "for project:", projects, "\n")
        log_message(warning_message, log_file)
      }
    }

    # Copy the relevant files
    files_for_phyloseq_object =
      list.files(source_folder,
                 pattern = "table.*\\.qza$|rooted-tree.*\\.qza$|classifier.*\\.qza|metadata\\.tsv$|metadata_extra\\.tsv$|dna-sequences.*\\.csv$|fcm.*\\.csv$|qPCR.*\\.csv$|prediction*\\.RDS$",
                 full.names = TRUE)

    file.copy(files_for_phyloseq_object, destination_folder, overwrite = TRUE)
  }
