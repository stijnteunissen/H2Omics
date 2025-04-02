barplot <- function(physeq = rarefied_genus_psmelt,
                    ntaxa = NULL,
                    norm_method = NULL,
                    sample_matrix = NULL,
                    group_by_factor = NULL,
                    taxrank = "Tax_label",
                    project_base_path = "/content/drive/MyDrive/H2Omics_workshop") {

  # Zorg dat benodigde libraries beschikbaar zijn
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  suppressMessages(library(ggplot2))

  # Convert copy_correction to lowercase for robust comparison
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder based on norm_method and copy_correction value
  destination_folder <- file.path(project_base_path, "sequencing_data", norm_method, cc_val)

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

  # Output folder for figures
  figures_folder <- file.path(project_base_path, "figures", "barplot")
  if (!dir.exists(figures_folder)) dir.create(figures_folder, recursive = TRUE)

  # Functie om PDF om te zetten naar PNG en weer te geven in Colab
  display_pdf_as_png <- function(pdf_file, output_png) {
    if (file.exists(pdf_file)) {
      pdf_convert_cmd <- sprintf("convert -density 150 %s -quality 90 %s", shQuote(pdf_file), shQuote(output_png))
      system(pdf_convert_cmd)

      if (file.exists(output_png)) {
        message("Converted to PNG: ", output_png)

        # Gebruik Python om afbeelding weer te geven in Colab
        py_run_string(sprintf("from IPython.display import display, Image; display(Image(filename='%s'))", output_png))
      } else {
        message("Failed to convert PDF to PNG: ", output_png)
      }
    } else {
      message("File does not exist: ", pdf_file)
    }
  }

  # List and process PDF files for relative plots
  pdf_files_relative <- list.files(barplot_folder, pattern = "barplot_relative\\.pdf$", full.names = TRUE)
  if (length(pdf_files_relative) > 0) {
    for (pdf_file in pdf_files_relative) {
      output_png <- file.path(figures_folder, paste0(basename(tools::file_path_sans_ext(pdf_file)), ".png"))
      display_pdf_as_png(pdf_file, output_png)
    }
  } else {
    message("No relative barplot PDFs found in ", barplot_folder)
  }

  # List and process PDF files for absolute plots
  pdf_files_absolute <- list.files(barplot_folder, pattern = "barplot_absolute\\.pdf$", full.names = TRUE)
  if (length(pdf_files_absolute) > 0) {
    for (pdf_file in pdf_files_absolute) {
      output_png <- file.path(figures_folder, paste0(basename(tools::file_path_sans_ext(pdf_file)), ".png"))
      display_pdf_as_png(pdf_file, output_png)
    }
  } else {
    message("No absolute barplot PDFs found in ", barplot_folder)
  }
}
