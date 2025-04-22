#' @export
psdata_to_tibble = function(physeq = rarefied_genus_physeq,
                            norm_method = NULL,
                            taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label")) {

  log_message(paste("Step 10: Convert phyloseq to tibble.", paste(projects, collapse = ", ")), log_file)

  message("phyloseq object transformed to a tibble")

  log_message("Successfully converted", log_file)
}
