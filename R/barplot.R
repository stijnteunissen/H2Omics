barplot = function(physeq = rarefied_genus_psmelt,
                   ntaxa = NULL,
                   norm_method = NULL,
                   sample_matrix = NULL,
                   group_by_factor = NULL,
                   taxrank = "Tax_label") {

  # Convert copy_correction to lowercase for robust comparison
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder based on norm_method and copy_correction value
  destination_folder <- file.path("/content/drive/MyDrive/H2Omics_workshop/sequencing_data", norm_method, cc_val)

  # Determine run option based on norm_method and present_variable_factors
  run_option <- NA
  if (norm_method == "fcm") {
    if (present_variable_factors == "treatment, timepoint") {
      run_option <- "run_option_1"
    } else if (present_variable_factors == "timepoint, treatment") {
      run_option <- "run_option_2"
    } else if (present_variable_factors == "timepoint") {
      run_option <- "run_option_3"
    } else if (present_variable_factors == "treatment") {
      run_option <- "run_option_4"
    }
  } else if (norm_method == "qpcr") {
    if (present_variable_factors == "sample_type, regrowth_day") {
      run_option <- "run_option_1"
    } else if (present_variable_factors == "regrowth_day, sample_type") {
      run_option <- "run_option_2"
    } else if (present_variable_factors == "sample_type") {
      run_option <- "run_option_3"
    } else if (present_variable_factors == "regrow_type") {
      run_option <- "run_option_4"
    }
  }

  if (is.na(run_option)) {
    stop("No valid run option determined based on the present_variable_factors.")
  }

  # Construct the folder path where barplot PDFs are stored
  barplot_folder <- file.path(destination_folder, run_option)
  message("Looking for barplot PDFs in: ", barplot_folder)

  # List PDF files in the barplot_folder
  pdf_files_relative <- list.files(barplot_folder, pattern = "barplot_relative\\.pdf$", full.names = TRUE)
  pdf_files_absolute <- list.files(barplot_folder, pattern = "barplot_absolute\\.pdf$", full.names = TRUE)

  display_pdf(pdf_files_relative)
  display_pdf(pdf_files_absolute)
}

