#' @export
beta_diversity <- function(physeq = physeq,
                           taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label"),
                           norm_method = NULL,
                           ordination_method = "PCoA",
                           color_factor = NULL,
                           color_continuous = FALSE,
                           shape_factor = NULL,
                           size_factor = NULL,
                           alpha_factor = NULL) {

  log_message(paste("Step 14: Making beta diversity.", paste(projects, collapse = ", ")), log_file)

  # Convert copy_correction to lowercase for robust comparison
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder based on norm_method and copy_correction value
  destination_folder <- file.path("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data", norm_method, cc_val, "After_cleaning_rds_files", taxrank_beta_div)

  psdata_file <- list.files(destination_folder, pattern = "corrected_counts\\.rds$", full.names = TRUE)

  psdata_relative = readRDS(psdata_file)

  figure_folder = paste0(base_path , projects, "/figures/Beta_diversity/")
  if (!dir.exists(figure_folder)) {
    dir.create(figure_folder, recursive = TRUE)
  }

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

  # Transform counts to relative abundance (percentage)
  psdata_relative <- transform_sample_counts(psdata_relative, function(x) x / sum(x) * 100)

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

  figure_file_path <- paste0(figure_folder, projects, "_beta_diversity_relative_", ordination_method, "_", taxrank_beta_div, "_level.pdf")
  ggsave(filename = figure_file_path, plot = combined_plot_relative, width = 10, height = 5)

  log_message("Beta diversity successfully plotted.", log_file)
}
