##########################
### PROCESSING DADA2 ####
########################

# List and write out sample names
get.sample.name <- function(fname) {
  base <- tools::file_path_sans_ext(basename(fname))
  sub("(_.*)?(-[^-]+)$", "", base)
}

### Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dnaTrnL <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orientsTrnL <- c(Forward = dnaTrnL, Complement = Biostrings::complement(dnaTrnL), Reverse = Biostrings::reverse(dnaTrnL), 
                   RevComp = Biostrings::reverseComplement(dnaTrnL))
  return(sapply(orientsTrnL, toString))  # Convert back to character vector
}

### Count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

### PRIMER OCCURENCE
primer_occurence <- function(fnFs, fnRs, FWD, REV){
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]])) 
}

### CUTADAPT
run_cutadapt <- function(i) {
  system2(
    cutadapt, args = c(
      R1.flags, R2.flags, 
      "-n", 2, 
      "-m", 21, '-M', 300, # see https://github.com/benjjneb/dada2/issues/2045#issuecomment-2449416862
      "-o", fnFs.cut[i], 
      "-p", fnRs.cut[i],
      fnFs.filtN[i], fnRs.filtN[i])
  )
}

### READS TRACKING
getN <- function(x) sum(getUniques(x))

transform_df <- function(object) { 
  require('dplyr')
  as.data.frame(object) %>% # transform an object/vector into a data.frame 
    rownames_to_column('Sample') # to perform the left_join afterwards
}

insert_middle <- function(original, insert, position) {
  require('tidyverse')
  str_c(str_sub(original, 1, position), insert, str_sub(original, position + 1, -1)) 
}

# original: The original string where we want to insert another string.
# insert: The string to inserted.
# position: The index (1-based) where the insert string will be added.


track_dada <- function(out.N, out,
                       dadaFs, dadaRs,
                       mergers,
                       seqtab.nochim) {
  require(dplyr, tibble, tidyr)
  
  out.N <- transform_df(out.N) # transform into data.frame
  out <- transform_df(out) %>% select(1,3) # selecting the 'Sample' and 'read.out' columns
  
  go.1 <- out.N %>%
    left_join(out, by = "Sample") %>%
    mutate(Sample = sub(insert_middle("-_.*$", barcode, 1), "", Sample)) %>% # enlever tout après ex: -ITS*
    mutate(Sample = paste0(year, '-', Sample)) # ajouter l'année
  
  # transform each variable
  go.2 <- transform_df(sapply(dadaFs, getN)) 
  go.3 <- transform_df(sapply(dadaRs, getN)) 
  go.4 <- transform_df(sapply(mergers_pooled, getN))
  go.5 <- transform_df(rowSums(seqtab.nochim))
  
  # combining each object 'go.X' to make a big table
  track <- left_join(go.1, go.2, by = "Sample") %>%
    left_join(., go.3, by = "Sample") %>% 
    left_join(., go.4, by = "Sample") %>% 
    left_join(., go.5, by = "Sample") %>% 
    replace(., is.na(.), 0)
  
  colnames(track) <- c("Sample", "input", "removeNs", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  
 # View result
  
  track %>% data.frame %>% 
    tibble %>% 
    mutate(N_filtering = (input-removeNs)/input,
           Quality_filtering = (removeNs-filtered)/removeNs,
           Denoising = (filtered-denoisedR)/filtered,
           Reads_merging = (denoisedR-merged)/denoisedR, # Proportion of reads lost to merging
           Bimera_removal = (merged-nonchim)/merged) %>% 
    pivot_longer(where(is.numeric), names_to = 'variable', values_to = 'values') %>% 
    mutate_all(~ ifelse(is.nan(.), 0, .)) # régler le problème des divsions par 0
}



###############
### PLOTS ######
#################

plot_track_change <- function(track_change) {
  require(dplyr, ggplot2)
  
  change_vars <- c('Bimera_removal', 'Reads_merging', 'Denoising', 'Quality_filtering', 'N_filtering')
  
  track_change %>% 
    filter(variable %in% change_vars) %>% 
    mutate(variable = factor(variable, level = change_vars)) %>% 
    ggplot(aes(y = variable, x = values)) +
    geom_boxplot() + theme_minimal() +
    labs(title = 'Proportion of reads lost at a specific pipeline step.')
}

topTaxa <- function(psmelt, taxLvl, topN) {
  psmelt %>% 
    group_by(!!sym(taxLvl)) %>% # group by tax level
    filter(relAb != 'NaN') %>% 
    summarise(relAb = mean(relAb)) %>% # find top abundant taxa
    arrange(desc(relAb)) %>% 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number() <= topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number() > topN ~ 'Others'))) # +1 to include the Others section!
}

