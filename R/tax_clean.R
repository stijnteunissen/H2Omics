#' @export
tax_clean = function(physeq = physeq, tax_filter = TRUE) {

  log_message(
    paste("Step 3: Tax clean: phyloseq taxa are cleaned.", paste(projects, collapse = ", ")), log_file)

  destination_folder <- file.path("/content/Workshop_H2Omics_test/H2Omics_workshop/sequencing_data/qpcr")

  uncleaned_file <- list.files(destination_folder, pattern = "phyloseq_uncleaned\\.rds$", full.names = TRUE)
  cleaned_file   <- list.files(destination_folder, pattern = "phyloseq_cleaned\\.rds$", full.names = TRUE)

  uncleaned_phyloseq <- readRDS(uncleaned_file)
  cleaned_phyloseq   <- readRDS(cleaned_file)

  keep_ranks <- c("Family", "Genus")

  subset_tax_table <- function(ps) {
    tax = tax_table(ps)
    ranks = intersect(colnames(tax), keep_ranks)
    tax_table(ps) = tax[, ranks, drop = FALSE]
    ps
  }

  uncleaned_phyloseq <- subset_tax_table(uncleaned_phyloseq)
  cleaned_phyloseq   <- subset_tax_table(cleaned_phyloseq)

  selected_otus <- c(
    "b79bafdef1d996f9c1c69b26b16a8d7d",
    "2b68abfcc52aa30e49215e5a7d0657dd",
    "bca1ec8c1c08341078a4a81852aec77c",
    "3203a668b6213d8fb17701f535e72334",
    "6a6110b7e3b01c8b2143982d71739c5b")

  uncleaned_sel <- prune_taxa(selected_otus, uncleaned_phyloseq)
  cleaned_sel   <- prune_taxa(selected_otus, cleaned_phyloseq)

  print(tax_table(uncleaned_sel))
  print(tax_table(cleaned_sel))

  new_folder    <- file.path(base_path, projects, "output_data/rds_files/Before_cleaning_rds_files")
  new_file_path <- file.path(new_folder, paste0(projects, "_phyloseq_cleaned.rds"))
  saveRDS(cleaned_phyloseq, file = new_file_path)

  message("Phyloseq object has been cleaned en opgeslagen.")
  log_message("Successfully Tax cleaned.", log_file)
}
