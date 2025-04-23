#' @export
heatmap <- function(physeq = rarefied_genus_psmelt,
                    ntaxa = NULL,
                    norm_method = NULL,
                    taxrank = "Tax_label") {

  log_message(paste("Step 12: making heatmap.", paste(projects, collapse = ", ")), log_file)

  # Convert copy_correction to lowercase for robust comparison
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder based on norm_method and copy_correction value
  destination_folder <- file.path("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data", norm_method, cc_val, "After_cleaning_rds_files/Tax_label")

  relative_files <- list.files(destination_folder, pattern = "pstibble_relative_data\\.rds$", full.names = TRUE)
  plot_data_rel = readRDS(relative_files)

  figure_folder = paste0(base_path , projects, "/figures/Heatmap/")
  if (!dir.exists(figure_folder)) {
    dir.create(figure_folder, recursive = TRUE)
  }

  base_heatmap = function(plot_data, x_value, abund_value, legend_name, x_label = "Sample") {
    ggplot(plot_data, aes(x = Sample, y = Tax_label)) +
      geom_tile(aes(fill = !!sym(abund_value)), color = NA) +
      scale_fill_gradient(low = "white", high = "darkred", name = legend_name) +
      labs(x = x_label, y = NULL) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_markdown(size = 10),
            axis.ticks.x = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            #strip.text = element_text(face = "bold", angle = 90, vjust = 0.5, hjust = 0),
            strip.text = element_text(face = "bold"),
            ggh4x.facet.nestline = element_line(colour = "black")) +
      scale_x_discrete(expand = c(0, 0)) +
      geom_text(aes(label = ifelse(mean_rel_abund > 3, paste0(round(mean_rel_abund, 0)), ifelse(mean_rel_abund == 0, ".", "")),
                    color = ifelse(mean_rel_abund > 50, "#D3D3D3", "black")),
                size = 3)
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

  variable_columns = intersect(present_variable_factors, colnames(plot_data_rel))
  factor_columns = unique(c(variable_columns))
  present_factors = if (length(factor_columns) > 0) factor_columns else NULL

  plot_data_rel_2 =
    plot_data_rel %>%
    group_by(Sample, Tax_label, na_type, !!!syms(present_factors)) %>%
    summarise(mean_rel_abund = sum(mean_rel_abund),
              .groups = "drop")

  if ("treatment" %in% colnames(plot_data_rel_2)) {
    plot_data_rel_2 <- plot_data_rel_2 %>%
      mutate(is_control = grepl("^untreated", tolower(treatment))) %>%
      mutate(Sample = factor(Sample, levels = unique(Sample[order(!is_control)])))
  }

  heatmap_relative =
    base_heatmap(plot_data_rel_2, "Sample", "mean_rel_abund", legend_name = "Relative\nAbundance (%)", x_label = "Sample") +
    scale_color_identity() +
    facet_add(present_factors)

  print(heatmap_relative)

  figure_file_path = paste0(figure_folder, projects, "_heatmap_relative.pdf")
  ggsave(filename = figure_file_path, plot = heatmap_relative, width = 12, height = 8, limitsize = FALSE)

  log_message("Heatmap successfully plotted.", log_file)
}
