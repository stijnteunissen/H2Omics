normalise_data <- function(physeq = without_mock_physeq,
                           norm_method = NULL,
                           copy_correction = TRUE) {

  # Convert copy_correction to a lowercase character string
  cc_val <- tolower(as.character(copy_correction))

  #destination_folder <- paste0("/wetsus_repo_analysis/r_visualisation_scripts/H2Omics_workshop/sequencing_data/", norm_method)
  destination_folder <- paste0("/content/drive/MyDrive/H2Omics_workshop/sequencing_data/", norm_method, "/", cc_val)

  if (cc_val == "false") {
    message("data normalised for copy number")
  } else if (cc_val == "false") {
    message("Only biomass normalisation applied for fcm.")
  }
}
