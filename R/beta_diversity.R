beta_diversity <- function(physeq = physeq,
                           taxrank = "ASV",
                           norm_method = NULL,
                           ordination_method = "PCoA",
                           color_factor = NULL,
                           color_continuous = TRUE,
                           shape_factor = NULL,
                           size_factor = NULL,
                           alpha_factor = NULL) {

  # Assume these globals are defined: copy_correction, base_path, projects, norm_method, taxrank
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder path
  # This folder should contain the pre-computed relative phyloseq object for the given taxrank
  destination_folder <- file.path("/content/drive/MyDrive/H2Omics_workshop/sequencing_data", norm_method, cc_val, taxrank)

  # Look for the relative phyloseq object in the destination folder
  ps_relative_files <- list.files(destination_folder, pattern = "physeq_relative\\.rds$", full.names = TRUE)
  if (length(ps_relative_files) == 0) {
    stop("No relative phyloseq object found in ", destination_folder)
  }
  psdata_relative <- readRDS(ps_relative_files[1])

  # Set up project and folder paths (these are used for saving figures)
  project_name <- projects
  project_folder <- file.path(base_path, project_name)
  figure_folder <- file.path(project_folder, "figures")
  beta_div_folder <- file.path(figure_folder, "Beta_diversity")
  if (!dir.exists(beta_div_folder)) { dir.create(beta_div_folder, recursive = TRUE) }

  # For a faster run time, if there are multiple na_types, subset to the first one.
  na_types <- unique(sample_data(psdata_relative)$na_type)
  if (length(na_types) > 1) {
    message("Multiple na_types detected; using the first: ", na_types[1])
    psdata_relative <- subset_samples(psdata_relative, na_type == na_types[1])
  }

  # Transform counts to relative abundance (if not already in percentage)
  psdata_relative <- transform_sample_counts(psdata_relative, function(x) x / sum(x) * 100)

  # Convert the specified variables to factors for plotting
  if (!is.null(color_factor)) {
    sample_data(psdata_relative)[[color_factor]] <- as.factor(sample_data(psdata_relative)[[color_factor]])
  }
  if (!is.null(shape_factor)) {
    sample_data(psdata_relative)[[shape_factor]] <- as.factor(sample_data(psdata_relative)[[shape_factor]])
  }
  if (!is.null(size_factor)) {
    sample_data(psdata_relative)[[size_factor]] <- as.factor(sample_data(psdata_relative)[[size_factor]])
  }
  if (!is.null(alpha_factor)) {
    sample_data(psdata_relative)[[alpha_factor]] <- as.factor(sample_data(psdata_relative)[[alpha_factor]])
  }

  # Define the base beta-diversity plotting function (this retains your formatting)
  base_beta_plot <- function(psdata, ordination_method, distance_method, title,
                             color_factor, shape_factor, size_factor, alpha_factor) {
    ordination_res <- ordinate(psdata, method = ordination_method, distance = distance_method)
    base_plot <- plot_ordination(psdata, ordination = ordination_res, axes = c(1,2))

    # Remove existing point layers to add custom ones
    if (length(base_plot$layers) > 0) {
      base_plot$layers <- base_plot$layers[!sapply(base_plot$layers, function(x) inherits(x$geom, "GeomPoint"))]
    }

    # Build aesthetic mapping
    aes_params <- list()
    if (!is.null(color_factor)) aes_params$color <- rlang::sym(color_factor)
    if (!is.null(shape_factor)) aes_params$shape <- rlang::sym(shape_factor)
    if (!is.null(size_factor))  aes_params$size  <- rlang::sym(size_factor)
    if (!is.null(alpha_factor)) aes_params$alpha <- rlang::sym(alpha_factor)

    base_plot <- base_plot + geom_point(mapping = do.call(aes, aes_params))

    # Apply color scales
    if (!is.null(color_factor)) {
      if (color_continuous) {
        base_plot <- base_plot + scale_color_continuous(low = "lightblue", high = "darkgreen")
      } else {
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
      base_plot <- base_plot + scale_alpha_continuous(range = c(0.3, 1))
    }

    base_plot <- base_plot +
      ggtitle(title) +
      labs(color = color_factor, shape = shape_factor, size = size_factor, alpha = alpha_factor) +
      theme(panel.background = element_rect(fill = "transparent"),
            panel.grid = element_line(colour = "grey90"),
            strip.text = element_text(face = "bold"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
            axis.line.y = element_line(color = "black", linewidth = 0.3),
            axis.line.x = element_line(color = "black", linewidth = 0.3)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

    return(base_plot)
  }

  # Create the ASV folder under beta-diversity figures
  asv_folder <- file.path(beta_div_folder, "ASV")
  if (!dir.exists(asv_folder)) { dir.create(asv_folder, recursive = TRUE) }

  # Generate beta diversity plots using different distance methods
  plot_Jac <- base_beta_plot(psdata_relative, ordination_method, "jaccard", "Jaccard (binary presence only)",
                             color_factor, shape_factor, size_factor, alpha_factor) + theme(legend.position = "none")
  plot_BC <- base_beta_plot(psdata_relative, ordination_method, "bray", "Bray-Curtis (presence + abundance)",
                            color_factor, shape_factor, size_factor, alpha_factor) + theme(legend.position = "none")
  plot_uu <- base_beta_plot(psdata_relative, ordination_method, "uunifrac", "Unweighted UniFrac (lineage presence only)",
                            color_factor, shape_factor, size_factor, alpha_factor) + theme(legend.position = "none")
  plot_wu <- base_beta_plot(psdata_relative, ordination_method, "wunifrac", "Weighted UniFrac (lineage presence and abundance)",
                            color_factor, shape_factor, size_factor, alpha_factor) + theme(legend.position = "none")

  # Extract legend and combine plots
  legend <- get_legend(plot_Jac + theme(legend.position = "right"))
  combined_plot_relative <- cowplot::plot_grid(plot_Jac, plot_BC, plot_uu, plot_wu, ncol = 2, labels = c("A", "B", "C", "D"))
  combined_plot_relative <- cowplot::plot_grid(combined_plot_relative, legend, ncol = 2, rel_widths = c(3, 0.8))

  print(combined_plot_relative)

  # Save the combined plot
  figure_file_path <- file.path(asv_folder, paste0(project_name, "_beta_diversity_relative_", ordination_method, "_", taxrank, "_level.pdf"))
  ggsave(filename = figure_file_path, plot = combined_plot_relative, width = 10, height = 5)
  log_message(paste("Relative beta diversity plot saved:", figure_file_path), log_file)

  return(combined_plot_relative)
}
