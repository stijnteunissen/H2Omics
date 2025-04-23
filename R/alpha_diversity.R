#' @export
alpha_diversity <- function(physeq = physeq,
                            norm_method = NULL,
                            taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label"),
                            date_factor = NULL) {

  log_message(paste("Step 13: Making alpha diversity.", paste(projects, collapse = ", ")), log_file)

  # Convert copy_correction to lowercase for robust comparison
  cc_val <- tolower(as.character(copy_correction))

  # Construct the destination folder based on norm_method and copy_correction value
  destination_folder <- file.path("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data", norm_method, cc_val, "After_cleaning_rds_files", "Tax_label")

  psdata_file <- list.files(destination_folder, pattern = paste0("phyloseq_Tax_label_level_", norm_method, "_normalised_cell_concentration_rarefied\\.rds$"), full.names = TRUE)

  psdata = readRDS(psdata_file)

  figure_folder = paste0(base_path , projects, "/figures/Alpha_diversity/")
  if (!dir.exists(figure_folder)) {
    dir.create(figure_folder, recursive = TRUE)
  }

  base_alpha_plot = function(alpha_data, x_value, y_value, x_label, y_label) {

    if (is.null(present_factors) || length(present_factors) == 0) {
      alpha_data <- alpha_data %>% mutate(grouping_factor = "all")
    } else {
      alpha_data <- alpha_data %>%
        mutate(grouping_factor = do.call(paste, c(across(all_of(present_factors)), sep = "_")))
    }

    plot = ggplot(alpha_data, aes(x = !!sym(x_value), y = !!sym(y_value), group = grouping_factor)) +
      #geom_jitter(aes(color = grouping_factor), size = 2, width = 0.2, show.legend = FALSE) +
      geom_col(fill = "steelblue", color = "steelblue", show.legend = FALSE) +
      #scale_color_manual(values = colorset) +
      theme_classic() +
      labs(x = x_label, y = y_label) +
      theme(
        legend.position = "bottom",
        legend.text = element_markdown(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        #strip.text = element_text(angle = 90, vjust = 0.5, hjust = 0),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black")) +
      scale_y_continuous(expand = c(0, 0))
  }

  facet_add = function(present_factors) {
    if (!is.null(present_factors) && length(present_factors) > 0) {
      return(facet_nested(cols = vars(!!!syms(present_factors)), scales = "free", space = "free", nest_line = element_line(linetype = 1)))
    } else {
      return(NULL)
    }
  }

  variable_columns = intersect(present_variable_factors, colnames(sample_data(psdata)))
  factor_columns = unique(c(variable_columns))
  present_factors = if (length(factor_columns) > 0) factor_columns else NULL

  alpha_data = suppressWarnings(estimate_richness(psdata, measures = c("Observed", "Chao1", "Shannon", "Simpson")))
  alpha_data = alpha_data %>% rownames_to_column(var = "sampleid")
  metadata = sample_data(psdata) %>% data.frame() %>% as_tibble()
  alpha_data_full = inner_join(metadata, alpha_data, by = "sampleid")

  if ("treatment" %in% colnames(alpha_data_full)) {
    alpha_data_full <- alpha_data_full %>%
      mutate(is_control = grepl("^untreated", treatment, ignore.case = TRUE)) %>%
      mutate(treatment = factor(treatment, levels = unique(treatment[order(!is_control)]))) %>%
      mutate(sampleid = factor(sampleid, levels = unique(sampleid[order(!is_control)])))
  }

  chao1_plot = base_alpha_plot(alpha_data_full, "sampleid", "Chao1", x_label = "Sample", y_label = "Chao1 Index") +
    facet_add(present_factors)

  shannon_plot = base_alpha_plot(alpha_data_full, "sampleid", "Shannon", x_label = "Sample", y_label = "Shannon Index") +
    facet_add(present_factors)

  combined_plot = plot_grid(chao1_plot + theme(legend.position = "none"), shannon_plot + theme(legend.position = "none"),
                            align = "hv", labels = c("A", "B"), nrow = 1)

  print(combined_plot)

  figure_file_path = paste0(figure_folder, projects, "_Genus_level_alpha_diversity.pdf")
  ggsave(filename = figure_file_path, plot = combined_plot, width = 12, height = 5)

  log_message("Alpha diversity successfully plotted.", log_file)
}
