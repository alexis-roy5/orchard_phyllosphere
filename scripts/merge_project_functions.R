####################################
######## MERGING TEST ##############
####################################

library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel)

# DESCRIPTION OF THE FUNCTION:
# seqtab building
  # Remove chimeras
# taxonomic assignement

# WRITE OUT sequence table & taxonomy table

# RETURN seqtab without chimeras

write.seqtab.taxonomic <- function(input.mergers, path.taxonomic, ref.db.path, ncores.used) { 
  # input.mergers: is the different options of merging
  # path.tax: where you want your .RDS
  # ref.db.path: path of your reference database
  # ncores : number of cores used by your computer
  
  # create the sequence table
  seqtab <- makeSequenceTable(input.mergers)
  
  # remove bimera
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = ncores.used, verbose = TRUE) 
  
  # look at the dimensions
  # print(dim(seqtab))
  # print(dim(seqtab.nochim))
  

  obj_name <- deparse(substitute(input.mergers))  # extract the name from input.mergers
  
  # adding in the middle the option NAME
  path.pre.seqRDS <- paste0(path.taxonomic, '/seqtab.')
  path.seqRDS <- paste0(path.pre.seqRDS, obj_name)
  
  write_rds(seqtab.nochim, paste0(path.seqRDS, '.RDS'))
  

  # creating the taxonomic table
  taxa <- assignTaxonomy(
    seqtab.nochim, 
    ref.db.path, 
    multithread= min(ncores.used, 72), tryRC = TRUE, verbose = TRUE)
  
  taxa_fixed <- taxa %>% as.data.frame() %>%
    mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .))) # modifying the table
  
  # adding in the middle the option NAME
  path.pre.taxRDS <- paste0(path.taxonomic, '/taxonomy.')
  path.taxRDS <- paste0(path.pre.taxRDS, obj_name)
  
  
  write_rds(taxa_fixed, paste0(path.taxRDS,'.RDS'))
  
  return(seqtab.nochim) # for further analysis
}

################################
######### PHYLOSEQ #############
################################

# PHYLOSEQ FUNCTIONS to run phyloseq.object

###############
# Functions ####
#################

# Keep samples with metadata
subset_samples <- function(seqtab, samples) {
  seqtab[rownames(seqtab) %in% samples, ] %>% # subset
    .[, colSums(.) > 0] # Remove ASVs with no hits
}

### %in% : est-ce que l'on trouve ce qu'il y a dans vecteur 1 dans vecteur 2 ?

# ASVs classified at the kingdom level and present in seqtab
subset_asvs <- function(taxonomy, seqtab, min_seq) {
  if (!is.data.frame(taxonomy)) {
    taxonomy <- as.data.frame(taxonomy)
  }
  
  asvs <- subset(taxonomy, Kingdom != "Unclassified") %>% # subset needs the input to be a df
    rownames %>% 
    intersect(
      colnames(seqtab)[colSums(seqtab) >= min_seq] 
    ) # only keep asvs still present in seqtab
  taxonomy[asvs, ] %>% as.matrix()
}

# Remove samples with fewer than n sequences once taxa removed

remove_ultra_rare <- function(seqtab, taxonomy, n) { 
  result <- seqtab[, rownames(taxonomy), drop = FALSE]  # Ensure it stays a data frame
  result <- result[rowSums(result) > n, , drop = FALSE]  # Filter rows (samples). n = sum across ASVs in a sample
  
  return(result)
}
###############################################

phyloseq.object <- function(seqtab, taxonomy, metadata, sample.names) {
  # seqtab: sequence table
  # taxonomy: taxonomic table
  # metadata: all the metadata that you want in your phyloseq object
  # sample.names : samples name
  require(phyloseq)
  
  # Keep samples with metadata info
  seqtab_ITS_sam <- subset_samples(seqtab, sample.names.phylo)
  
  # Subset ASVs (ex: 100 hits)
  taxa_ITS_sam <- subset_asvs(taxonomy, seqtab_ITS_sam, 100)
  
  # Remove near-empty samples - LOOK at the number and change it in consequences
  seqtab_ITS_sam_filt <- remove_ultra_rare(seqtab_ITS_sam, taxa_ITS_sam, 3000) 
  
  # Phyloseq objects
  ps_ITS <- phyloseq(
    tax_table(taxa_ITS_sam),
    otu_table(seqtab_ITS_sam_filt, taxa_are_rows = FALSE),
    sample_data(samples)
  )
  
  saveRDS(ps_ITS, paste0(path.out.ps, "ps_ITS.rds"))
}





