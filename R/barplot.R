barplot <- function(physeq = rarefied_genus_psmelt,
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

  # Function to display a PDF file using IRdisplay by reading its raw content
  display_pdf_file <- function(pdf_file) {
    if (file.exists(pdf_file)) {
      pdf_size <- file.info(pdf_file)$size
      pdf_data <- readBin(pdf_file, what = "raw", n = pdf_size)
      IRdisplay::display_pdf(pdf_data)
    } else {
      message("File does not exist: ", pdf_file)
    }
  }

  # List and display PDF files for relative plots
  pdf_files_relative <- list.files(barplot_folder, pattern = "barplot_relative\\.pdf$", full.names = TRUE)
  if (length(pdf_files_relative) > 0) {
    for (pdf_file in pdf_files_relative) {
      display_pdf_file(pdf_file)
    }
  } else {
    message("No relative barplot PDFs found in ", barplot_folder)
  }

  # List and display PDF files for absolute plots
  pdf_files_absolute <- list.files(barplot_folder, pattern = "barplot_absolute\\.pdf$", full.names = TRUE)
  if (length(pdf_files_absolute) > 0) {
    for (pdf_file in pdf_files_absolute) {
      display_pdf_file(pdf_file)
    }
  } else {
    message("No absolute barplot PDFs found in ", barplot_folder)
  }
}
