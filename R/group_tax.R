#' @export
group_tax = function(physeq = rarefied_asv_physeq,
                     norm_method = NULL,
                     taxrank = c("Phylum", "Class", "Order", "Family", "Tax_label")) {
  log_message(paste("Step 9: Tax glom: OTUs are merged at different taxonomic levels.", paste(projects, collapse = ", ")), log_file)

  message("features grouped on taxonomic ranks")

  log_message("Successfully Tax glom.", log_file)
}
