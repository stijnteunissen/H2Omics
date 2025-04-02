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

  relative_files <- list.files(destination_folder, pattern = "relative_data\\.rds$", full.names = TRUE)
  absolute_files <- list.files(destination_folder, pattern = "absolute_data\\.rds$", full.names = TRUE)

  plot_data_rel = readRDS(relative_files)
  plot_data_norm = readRDS(absolute_files)

  # Helper function to generate the base barplot
  base_barplot <- function(plot_data, x_value, y_value, colorset,
                           x_label = "Sample", y_label = "Value") {
    p <- ggplot(plot_data, aes(x = !!sym(x_value), y = !!sym(y_value), fill = Tax_label)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name = "Genus", values = colorset) +
      theme_classic(base_size = 14) +
      labs(x = x_label, y = y_label, fill = "Genus") +
      theme(axis.ticks.x = element_blank(),
            legend.text = element_markdown(),
            legend.key.size = unit(5, "pt"),
            legend.position = "bottom",
            strip.background = element_rect(colour = "white"),
            strip.text = element_text(face = "bold"),
            ggh4x.facet.nestline = element_line(colour = "black"))

    if (!is.null(present_factors)) {
      p <- p + theme(axis.text.x = element_blank())
    } else {
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
    }
    return(p)
  }

  # Helper function to add facets if present_factors is provided
  facet_add <- function(present_factors) {
    if (!is.null(present_factors) && length(present_factors) > 0) {
      return(
        facet_nested(
          cols = vars(!!!syms(present_factors)),
          scales = "free_x",
          space = "free",
          nest_line = element_line(linetype = 1)
        )
      )
    } else {
      return(NULL)
    }
  }

  # Set default color palette if not provided
  if (is.null(colorset)) {
    dark2_colors <- brewer.pal(8, "Dark2")
    paired_colors <- brewer.pal(12, "Paired")
    set_colors <- brewer.pal(8, "Set2")
    set1_colors <- brewer.pal(8, "Set1")
    spectral_colors <- brewer.pal(11, "Spectral")
    additional_palette <- unique(c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf",
      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
      "#393b79", "#9c755f", "#e7298a", "#66c2a5", "#fc8d62",
      "#8da0cb", "#e78ac3", "#a6d854", "#ffed6f", "#ffff33", "#fdbf6f", "#ff7f00", "#6a3d9a", "#b15928",
      "#e41a1c", "#377eb8", "#4daf4a", "#ff6a4d",
      "#c6dbef", "#fdae61", "#fee08b", "#91bfdb", "#d73027", "#4575b4", "#313695", "#ffcc00", "#a1d99b",
      "#ff99cc", "#32cd32", "#ff6347", "#20b2aa", "#c71585", "#3cb371", "#6495ed", "#9b59b6",
      "#2ecc71", "#e74c3c", "#3498db", "#f39c12", "#8e44ad", "#16a085", "#f1c40f", "#d35400", "#27ae60",
      "#2980b9", "#e67e22"
    ))
    colorset <- unique(c(dark2_colors, paired_colors, set_colors, set1_colors, spectral_colors, additional_palette))
  }

  # Set y-axis labels for relative and absolute plots
  ylabel_rel <- "Relative Abundance (%)"
  ylabel_abs <- ifelse(sample_matrix == "liquid",
                       "Cell equivalents (Cells/ml) sample",
                       "Cell equivalents (Cells/gram) sample")

  # Create output folder if it does not exist
  barplot_folder <- file.path(figure_folder, "Barplot")
  if (!dir.exists(barplot_folder)) {
    dir.create(barplot_folder, recursive = TRUE)
  }

  ### Generate Relative Plot ###
  if (!is.null(group_by_factor)) {
    all_plots_rel <- list()
    factors <- unique(plot_data_rel[[group_by_factor]])
    for (x in factors) {
      data_filtered <- plot_data_rel %>% filter(.data[[group_by_factor]] == x)
      plot_rel <- base_barplot(data_filtered, "Sample", "mean_rel_abund", colorset,
                               x_label = "Sample", y_label = ylabel_rel) +
        facet_add(present_factors) +
        scale_y_continuous(expand = c(0, 0))
      all_plots_rel[[x]] <- plot_rel
    }
    legend_rel <- get_legend(all_plots_rel[[1]] + theme(legend.position = "right"))
    nplots_rel <- length(all_plots_rel)
    final_plot_rel <- wrap_plots(all_plots_rel, ncol = nplots_rel) +
      plot_annotation(title = "Relative Abundance")
    barplot_relative <- plot_grid(final_plot_rel, legend_rel, ncol = 1, rel_heights = c(3, 1))
  } else {
    barplot_relative <- base_barplot(plot_data_rel, "Sample", "mean_rel_abund", colorset,
                                     x_label = "Sample", y_label = ylabel_rel) +
      facet_add(present_factors) +
      scale_y_continuous(expand = c(0, 0))
  }

  print(barplot_relative)
  rel_file_pdf <- file.path(barplot_folder, paste0(project_name, "_barplot_relative.pdf"))
  ggsave(filename = rel_file_pdf, plot = barplot_relative, width = 12, height = 8)

  ### Generate Absolute Plot ###
  if (!is.null(group_by_factor)) {
    all_plots_abs <- list()
    factors <- unique(plot_data_abs[[group_by_factor]])
    for (x in factors) {
      data_filtered <- plot_data_abs %>% filter(.data[[group_by_factor]] == x)
      plot_abs <- base_barplot(data_filtered, "Sample", "norm_abund", colorset,
                               x_label = "Sample", y_label = ylabel_abs) +
        facet_add(present_factors) +
        scale_y_continuous(expand = c(0, 0))
      all_plots_abs[[x]] <- plot_abs
    }
    legend_abs <- get_legend(all_plots_abs[[1]] + theme(legend.position = "right"))
    nplots_abs <- length(all_plots_abs)
    final_plot_abs <- wrap_plots(all_plots_abs, ncol = nplots_abs) +
      plot_annotation(title = "Absolute Abundance")
    barplot_absolute <- plot_grid(final_plot_abs, legend_abs, ncol = 1, rel_heights = c(3, 1))
  } else {
    barplot_absolute <- base_barplot(plot_data_abs, "Sample", "norm_abund", colorset,
                                     x_label = "Sample", y_label = ylabel_abs) +
      facet_add(present_factors) +
      scale_y_continuous(expand = c(0, 0))
  }

  print(barplot_absolute)
  abs_file_pdf <- file.path(barplot_folder, paste0(project_name, "_barplot_absolute.pdf"))
  ggsave(filename = abs_file_pdf, plot = barplot_absolute, width = 12, height = 8)
}
