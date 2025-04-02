beta_diversity <- function(physeq = physeq,
                           taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label"),
                           norm_method = NULL,
                           ordination_method = "PCoA",
                           color_factor = NULL,
                           color_continuous = TRUE,
                           shape_factor = NULL,
                           size_factor = NULL,
                           alpha_factor = NULL) {

  # Set up project and folder paths
  project_name <- projects
  project_folder <- paste0(base_path, project_name)
  figure_folder <- paste0(project_folder, "/figures/")
  destination_folder <- paste0(project_folder, "/input_data/")
  output_folder_csv_files <- paste0(project_folder, "/output_data/csv_files/")

  # Function for creating the base beta-diversity plot
  base_beta_plot <- function(psdata, ordination_method, distance_method, title,
                             color_factor, shape_factor, size_factor, alpha_factor) {
    # Perform the ordination
    ordination_res <- ordinate(psdata, method = ordination_method, distance = distance_method)
    base_plot <- plot_ordination(psdata, ordination = ordination_res, axes = c(1, 2))

    # Remove any existing point layers to add our own
    if (length(base_plot$layers) > 0) {
      base_plot$layers <- base_plot$layers[!sapply(base_plot$layers, function(x) inherits(x$geom, "GeomPoint"))]
    }

    # Build the aesthetic mapping list
    aes_params <- list()
    if (!is.null(color_factor)) aes_params$color <- sym(color_factor)
    if (!is.null(shape_factor)) aes_params$shape <- sym(shape_factor)
    if (!is.null(size_factor))  aes_params$size  <- sym(size_factor)
    if (!is.null(alpha_factor)) aes_params$alpha <- sym(alpha_factor)

    # Add points with the specified aesthetics
    base_plot <- base_plot + geom_point(mapping = do.call(aes, aes_params))

    # Apply color scales based on whether the color factor is continuous or discrete
    if (!is.null(color_factor)) {
      if (color_continuous == TRUE) {
        base_plot <- base_plot + scale_color_continuous(low = "lightblue", high = "darkgreen")
      } else if (color_continuous == FALSE) {
        base_plot <- base_plot + scale_color_manual(values = colorset)
      }
    }
    if (!is.null(shape_factor)) {
      base_plot <- base_plot + scale_shape_manual(values = shapeset)
    }
    if (!is.null(size_factor)) {
      base_plot <- base_plot + scale_size_manual(values = sizeset)
    }
    if (!is.null(alpha_factor)) {
      base_plot <- base_plot + scale_alpha_continuous(range = c(0.3, 1))  # Set transparency for alpha
    }

    # Add title and labels; adjust theme and axis settings
    base_plot <- base_plot +
      ggtitle(title) +
      labs(color = color_factor, shape = shape_factor, size = size_factor, alpha = alpha_factor) +
      theme(panel.background = element_rect(fill = "transparent"),
            panel.grid = element_line(colour = "grey90"),
            strip.text = element_text(face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
            axis.line.y = element_line(color = "black", linewidth = 0.3),
            axis.line.x = element_line(color = "black", linewidth = 0.3)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))  # 5% expansion for y-axis

    return(base_plot)
  }

  # Convert copy_correction to lowercase for robust comparison
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder based on norm_method and copy_correction value
  destination_folder <- file.path("/content/drive/MyDrive/H2Omics_workshop/sequencing_data", norm_method, cc_val)

  # Create the main beta-diversity folder
  beta_div_folder <- paste0(figure_folder, "Beta_diversity/")
  if (!dir.exists(beta_div_folder)){dir.create(beta_div_folder)}

  if (tolower(taxrank[1]) == "asv") {
    log_message("Processing ASV-level beta diversity", log_file)

    if (norm_method == "fcm") {
      psdata_relative_file <- readRDS()
    } else if (norm_method == "qpcr") {
      psdata_relative_file <- physeq[["psdata_asv_copy_number_corrected"]]
    }

    psdata_relative = readRDS(psdata_relative_file)

    # Create the ASV folder under beta-diversity
    asv_folder <- paste0(beta_div_folder, "ASV/")
    if (!dir.exists(asv_folder)){dir.create(asv_folder)}

    # Transform counts to relative abundance (percentage)
    psdata_relative <- transform_sample_counts(psdata_relative, function(x) x / sum(x) * 100)

    # Convert specified variables to factors for absolute data (if available)
    if (!is.null(norm_method)) {
      if (!is.null(color_factor))
        sample_data(psdata_absolute)[[color_factor]] <- as.factor(sample_data(psdata_absolute)[[color_factor]])
      if (!is.null(shape_factor))
        sample_data(psdata_absolute)[[shape_factor]] <- as.factor(sample_data(psdata_absolute)[[shape_factor]])
      if (!is.null(size_factor))
        sample_data(psdata_absolute)[[size_factor]] <- as.factor(sample_data(psdata_absolute)[[size_factor]])
      if (!is.null(alpha_factor))
        sample_data(psdata_absolute)[[alpha_factor]] <- as.factor(sample_data(psdata_absolute)[[alpha_factor]])
    }

    # Define color, shape, and size sets based on psdata_relative
    if (!is.null(color_factor)) {
      sample_data(psdata_relative)[[color_factor]] <- as.factor(sample_data(psdata_relative)[[color_factor]])
      unique_colors <- levels(sample_data(psdata_relative)[[color_factor]])
      n_colors <- length(unique_colors)
      if (color_continuous == FALSE) {
        colorset <<- scales::hue_pal()(n_colors)
      }
    }
    if (!is.null(shape_factor)) {
      sample_data(psdata_relative)[[shape_factor]] <- as.factor(sample_data(psdata_relative)[[shape_factor]])
      n_shapes <- length(unique(sample_data(psdata_relative)[[shape_factor]]))
      shapeset <<- seq_len(n_shapes)
    }
    if (!is.null(size_factor)) {
      sample_data(psdata_relative)[[size_factor]] <- as.factor(sample_data(psdata_relative)[[size_factor]])
      n_sizes <- length(unique(sample_data(psdata_relative)[[size_factor]]))
      sizeset <<- seq(2, 2 + 1.2 * (n_sizes - 1), by = 1.2)
    }
    if (!is.null(alpha_factor)) {
      sample_data(psdata_relative)[[alpha_factor]] <- as.factor(sample_data(psdata_relative)[[alpha_factor]])
      alphaset <<- scale_alpha_continuous(range = c(0.3, 1))
    }

      # Create beta-diversity plots for a single na_type
      plot_Jac <- base_beta_plot(psdata_relative, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                 color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_BC <- base_beta_plot(psdata_relative, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_uu <- base_beta_plot(psdata_relative, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")
      plot_wu <- base_beta_plot(psdata_relative, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                color_factor, shape_factor, size_factor, alpha_factor) +
        theme(legend.position = "none")

      legend <- get_legend(plot_Jac + theme(legend.position = "right"))
      combined_plot_relative <- cowplot::plot_grid(plot_Jac, plot_BC, plot_uu, plot_wu,
                                                   ncol = 2, labels = c("A", "B", "C", "D"))
      combined_plot_relative <- cowplot::plot_grid(combined_plot_relative, legend, ncol = 2,
                                                   rel_widths = c(3, 0.8))

      print(combined_plot_relative)

      # Save the relative beta diversity plot (using the provided 'level' in the filename)
      figure_file_path <- paste0(asv_folder, project_name, "_beta_diversity_relative_", ordination_method, "_asv_level.pdf")
      ggsave(filename = figure_file_path, plot = combined_plot_relative, width = 10, height = 5)
      log_message(paste("Relative beta diversity plot saved:", figure_file_path), log_file)

  } else {
    for (tax in taxrank) {
      log_message(paste("Processing taxonomic level:", tax), log_file)

      if (is.null(norm_method)) {
        psdata_relative <- physeq
        psdata_absolute <- NULL
      } else if (norm_method == "fcm") {
        psdata_relative <- physeq[[paste0("psdata_copy_number_corrected_", tax)]]
        psdata_absolute <- physeq[[paste0("psdata_fcm_norm_rarefied_", tax)]]
      } else if (norm_method == "qpcr") {
        psdata_relative <- physeq[[paste0("psdata_copy_number_corrected_", tax)]]
        psdata_absolute <- physeq[[paste0("psdata_qpcr_norm_rarefied_", tax)]]
      }

      tax_folder <- paste0(beta_div_folder, tax, "/")
      if (!dir.exists(tax_folder)) {
        dir.create(tax_folder, recursive = TRUE)
      }

      psdata_relative <- transform_sample_counts(psdata_relative, function(x) x / sum(x) * 100)

      if (!is.null(norm_method)) {
        if (!is.null(color_factor))
          sample_data(psdata_absolute)[[color_factor]] <- as.factor(sample_data(psdata_absolute)[[color_factor]])
        if (!is.null(shape_factor))
          sample_data(psdata_absolute)[[shape_factor]] <- as.factor(sample_data(psdata_absolute)[[shape_factor]])
        if (!is.null(size_factor))
          sample_data(psdata_absolute)[[size_factor]] <- as.factor(sample_data(psdata_absolute)[[size_factor]])
        if (!is.null(alpha_factor))
          sample_data(psdata_absolute)[[alpha_factor]] <- as.factor(sample_data(psdata_absolute)[[alpha_factor]])
      }

      if (!is.null(color_factor)) {
        sample_data(psdata_relative)[[color_factor]] <- as.factor(sample_data(psdata_relative)[[color_factor]])
        unique_colors <- levels(sample_data(psdata_relative)[[color_factor]])
        n_colors <- length(unique_colors)
        if (color_continuous == FALSE) {
          colorset <<- scales::hue_pal()(n_colors)
        }
      }
      if (!is.null(shape_factor)) {
        sample_data(psdata_relative)[[shape_factor]] <- as.factor(sample_data(psdata_relative)[[shape_factor]])
        n_shapes <- length(unique(sample_data(psdata_relative)[[shape_factor]]))
        shapeset <<- seq_len(n_shapes)
      }
      if (!is.null(size_factor)) {
        sample_data(psdata_relative)[[size_factor]] <- as.factor(sample_data(psdata_relative)[[size_factor]])
        n_sizes <- length(unique(sample_data(psdata_relative)[[size_factor]]))
        sizeset <<- seq(2, 2 + 1.2 * (n_sizes - 1), by = 1.2)
      }
      if (!is.null(alpha_factor)) {
        sample_data(psdata_relative)[[alpha_factor]] <- as.factor(sample_data(psdata_relative)[[alpha_factor]])
        alphaset <<- scale_alpha_continuous(range = c(0.3, 1))
      }

        plot_Jac <- base_beta_plot(psdata_relative, ordination_method, "jaccard", "Jaccard\n(binary presence only)",
                                   color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_BC <- base_beta_plot(psdata_relative, ordination_method, "bray", "Bray-Curtis\n(presence + abundance)",
                                  color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_uu <- base_beta_plot(psdata_relative, ordination_method, "uunifrac", "Unweighted UniFrac\n(lineage presence only)",
                                  color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")
        plot_wu <- base_beta_plot(psdata_relative, ordination_method, "wunifrac", "Weighted UniFrac\n(lineage presence and abundance)",
                                  color_factor, shape_factor, size_factor, alpha_factor) +
          theme(legend.position = "none")

        legend <- get_legend(plot_Jac + theme(legend.position = "right"))
        combined_plot_relative <- cowplot::plot_grid(plot_Jac, plot_BC, plot_uu, plot_wu,
                                                     ncol = 2, labels = c("A", "B", "C", "D"))
        combined_plot_relative <- cowplot::plot_grid(combined_plot_relative, legend, ncol = 2,
                                                     rel_widths = c(3, 0.8))

        print(combined_plot_relative)

        figure_file_path <- paste0(tax_folder, project_name, "_beta_diversity_relative_", ordination_method, "_", tax, "_level.pdf")
        ggsave(filename = figure_file_path, plot = combined_plot_relative, width = 10, height = 5)
        log_message(paste("Relative beta diversity plot saved:", figure_file_path), log_file)
    }
  }
}



