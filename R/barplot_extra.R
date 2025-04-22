#' @export
barplot_extra = function(physeq = rarefied_genus_psmelt,
                   ntaxa = NULL,
                   norm_method = NULL,
                   sample_matrix = NULL,
                   group_by_factor = NULL,
                   taxrank = "Tax_label") {

  if (norm_method == "qpcr") {
  # Convert copy_correction to lowercase for robust comparison
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder based on norm_method and copy_correction value
  destination_folder <- file.path("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data", norm_method, cc_val, "After_cleaning_rds_files")

  relative_files <- list.files(destination_folder, pattern = "rel_plot_data\\.rds$", full.names = TRUE)
  absolute_files <- list.files(destination_folder, pattern = "norm_plot_data\\.rds$", full.names = TRUE)

  plot_data_rel = readRDS(relative_files)
  plot_data_norm = readRDS(absolute_files)

  figure_folder = paste0(base_path , projects, "/figures/")
  if (!dir.exists(figure_folder)) {
    dir.create(figure_folder, recursive = TRUE)
  }

  # base barplot
  base_barplot = function(plot_data, x_value, y_value, colorset, x_label = "Sample", y_label = "Cell equivalents (Cells/ml) sample") {
    p = ggplot(plot_data, aes(x = !!sym(x_value), y = !!sym(y_value), fill = Tax_label)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(name = "Genus", values = colorset) +
      theme_classic(base_size = 13) +
      labs(x = x_label, y = y_label, fill = "Genus") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.text = element_markdown(),
            legend.key.size = unit(5, "pt"),
            legend.position = "bottom",
            strip.background = element_rect(colour = "white"),
            strip.text = element_text(face = "bold"),
            #strip.text = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 0),
            ggh4x.facet.nestline = element_line(colour = "black")) +
      guides(fill = guide_legend(nrow = 8))
  }

  facet_add = function(present_factors) {
    if (!is.null(present_factors) && length(present_factors) > 0) {
      return(
        facet_nested(
          cols = vars(!!!syms(present_factors)), scales = "free_x", space = "free", nest_line = element_line(linetype = 1)))
    } else {
      return(NULL)
    }
  }

  # colorset
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

  variable_columns = intersect(present_variable_factors, colnames(plot_data_rel))
  factor_columns = unique(c(variable_columns))
  present_factors = if (length(factor_columns) > 0) factor_columns else NULL

  plot_data_rel %>%
    arrange(desc(mean_rel_abund)) %>%
    pull(Tax_label) %>%
    unique() -> names(colorset)

  if ("treatment" %in% colnames(plot_data_rel)) {
    plot_data_rel <- plot_data_rel %>%
      mutate(is_control = grepl("^untreated", tolower(treatment))) %>%
      mutate(Sample = factor(Sample, levels = unique(Sample[order(!is_control)])))
  }

  barplot_relative =
    base_barplot(plot_data_rel, "Sample", "mean_rel_abund", colorset, x_label = "Sample", y_label = "Relative Abundance (%)") +
    facet_add(present_factors) +
    scale_y_continuous(expand = c(0, 0))

  print(barplot_relative)

  figure_file_path = paste0(figure_folder, projects, "_barplot_pathogens_rel.pdf")
  ggsave(filename = figure_file_path, plot = barplot_relative, width = 12, height = 8)

  } else if (norm_method == "fcm") {
  message("This plot is only available for the other dataset. Return to chapter 2.1, start a new project, and re-run everything below.")
  }
}
