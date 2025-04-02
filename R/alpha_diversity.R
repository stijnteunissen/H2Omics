alpha_diversity <- function(physeq = physeq,
                            norm_method = NULL,
                            taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label"),
                            date_factor = NULL) {
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

  # Construct the folder path where alpha diversity PDFs are stored
  alpha_div_folder <- file.path(destination_folder, run_option, "Alpha_diversity", taxrank)
  message("Looking for alpha diversity PDFs in: ", alpha_div_folder)

  # List PDF files in the alpha_div_folder
  pdf_files <- list.files(alpha_div_folder, pattern = "alpha_diversity\\.pdf$", full.names = TRUE)

  if (length(pdf_files) == 0) {
    message("No alpha diversity PDFs found in ", alpha_div_folder)
  } else {
    # Install IRdisplay if not available
    if (!requireNamespace("IRdisplay", quietly = TRUE)) {
      install.packages("IRdisplay")
    }
    suppressMessages(library(IRdisplay))

    # Display each PDF inline
    for (pdf_file in pdf_files) {
      message("Displaying: ", pdf_file)
      pdf_size <- file.info(pdf_file)$size
      pdf_data <- readBin(pdf_file, what = "raw", n = pdf_size)
      IRdisplay::display_pdf(pdf_data)
    }
  }
}
