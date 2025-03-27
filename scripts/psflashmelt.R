# the original psmelt function, but optimized using data.table and dplyr syntax
# should be insanely faster than the original
# shamelessly AI generated
# Creates a much smaller table because it omits every row for which Abundance == 0

psmelt_dt <- function(physeq) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' needed for this function. Please install it.")
  }
  
  if (!inherits(physeq, "phyloseq")) {
    rankNames = NULL
    sampleVars = NULL
  }
  else {
    rankNames = rank_names(physeq, FALSE)
    sampleVars = sample_variables(physeq, FALSE)
  }
  reservedVarnames = c("Sample", "Abundance", "OTU")
  type1aconflict = intersect(reservedVarnames, sampleVars)
  if (length(type1aconflict) > 0) {
    wh1a = which(sampleVars %in% type1aconflict)
    new1a = paste0("sample_", sampleVars[wh1a])
    warning("The sample variables: \n", paste(sampleVars[wh1a], 
                                              collapse = ", "), "\n have been renamed to: \n", 
            paste0(new1a, collapse = ", "), "\n", "to avoid conflicts with special phyloseq plot attribute names.")
    colnames(sample_data(physeq))[wh1a] <- new1a
  }
  type1bconflict = intersect(reservedVarnames, rankNames)
  if (length(type1bconflict) > 0) {
    wh1b = which(rankNames %in% type1bconflict)
    new1b = paste0("taxa_", rankNames[wh1b])
    warning("The rank names: \n", paste(rankNames[wh1b], 
                                        collapse = ", "), "\n have been renamed to: \n", 
            paste0(new1b, collapse = ", "), "\n", "to avoid conflicts with special phyloseq plot attribute names.")
    colnames(tax_table(physeq))[wh1b] <- new1b
  }
  type2conflict = intersect(sampleVars, rankNames)
  if (length(type2conflict) > 0) {
    wh2 = which(sampleVars %in% type2conflict)
    new2 = paste0("sample_", sampleVars[wh2])
    warning("The sample variables: \n", paste0(sampleVars[wh2], 
                                               collapse = ", "), "\n have been renamed to: \n", 
            paste0(new2, collapse = ", "), "\n", "to avoid conflicts with taxonomic rank names.")
    colnames(sample_data(physeq))[wh2] <- new2
  }
  
  # 
  otutab = otu_table(physeq)
  if (!taxa_are_rows(otutab)) {
    otutab <- t(otutab)
  }
  
  # Convert to long format using tidyverse approach
  mdf <- as(otutab, "matrix") %>%
    tibble::as_tibble(rownames = "OTU") %>%
    tidyr::pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    dplyr::filter(Abundance > 0)
  
  # Sample data merge
  if (!is.null(sampleVars)) {
    sdf <- sample_data(physeq) %>% 
      as("data.frame") %>% 
      tibble::rownames_to_column("Sample")
    mdf <- dplyr::left_join(mdf, sdf, by = "Sample")
  }
  
  # Tax table merge
  if (!is.null(rankNames)) {
    TT = access(physeq, "tax_table")
    keepTTcols <- colSums(is.na(TT)) < ntaxa(TT)
    if (length(which(keepTTcols)) > 0 & ncol(TT) > 0) {
      tdf <- TT[, keepTTcols] %>%
        as("matrix") %>%
        tibble::as_tibble(rownames = "OTU")
      mdf <- dplyr::left_join(mdf, tdf, by = "OTU")
    }
  }
  
  mdf <- dplyr::arrange(mdf, desc(Abundance))
  return(as.data.frame(mdf))
}
