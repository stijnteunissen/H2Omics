normalise_data <- function(physeq = without_mock_physeq,
                           norm_method = NULL,
                           copy_correction = TRUE) {

  install.packages("IRdisplay", quiet = TRUE)
  suppressMessages(library(IRdisplay))

  # Convert copy_correction to a lowercase character string
  cc_val <- tolower(as.character(copy_correction))

  #destination_folder <- paste0("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method)
  destination_folder <- paste0("/content/drive/MyDrive/H2Omics_workshop/sequencing_data/", norm_method, "/", cc_val)

  if (cc_val == "false") {

    # List PDF files in the destination folder matching the pattern
    pdf_files <- list.files(destination_folder, pattern = "copy_number_comparison.*\\.pdf$", full.names = TRUE)

    pdf_file <- pdf_files

    IRdisplay::display_pdf(pdf_file)

  } else if (cc_val == "false") {
    message("Only biomass normalisation applied for fcm.")
  } else {
    message("data normalised for copy number")
  }
}
