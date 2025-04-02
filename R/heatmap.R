heatmap <- function(physeq = rarefied_genus_psmelt,
                    ntaxa = NULL,
                    norm_method = NULL,
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

  # Construct the folder path where heatmap PDFs are stored
  heatmap_folder <- file.path(destination_folder, run_option)
  message("Looking for heatmap PDFs in: ", heatmap_folder)

  # List PDF files in the heatmap_folder matching the pattern "heatmap_relative.pdf"
  pdf_files <- list.files(heatmap_folder, pattern = "heatmap_relative\\.pdf$", full.names = TRUE)

  # Display the first PDF file if available using IRdisplay::display_pdf
  if (length(pdf_files) == 0) {
    message("No heatmap PDFs found in ", heatmap_folder)
  } else {
    # Use the first matching file
    pdf_file <- pdf_files[1]
    pdf_size <- file.info(pdf_file)$size
    pdf_data <- readBin(pdf_file, what = "raw", n = pdf_size)
    if (!requireNamespace("IRdisplay", quietly = TRUE)) {
      install.packages("IRdisplay")
    }
    suppressMessages(library(IRdisplay))
    IRdisplay::display_pdf(pdf_data)
  }
}
