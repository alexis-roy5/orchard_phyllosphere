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

track_dada <- function(out.N, out,
                       dadaFs, dadaRs,
                       mergers,
                       seqtab.nochim) {
  
  track <- cbind(out.N, out[,2], 
                 sapply(dadaFs, getN), 
                 sapply(dadaRs, getN), 
                 sapply(mergers, getN),
                 rowSums(seqtab.nochim)) 
  # trouver le nombre de read unique -> ASVs
  
  colnames(track) <- c("input", "removeNs", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  
  track %>% data.frame %>% 
    rownames_to_column('Sample') %>% 
    tibble %>% 
    filter(filtered>10) %>% 
    mutate(lost_Ns = (input-removeNs)/input,
           lost_filt = (removeNs-filtered)/removeNs,
           lost_noise = (filtered-denoisedR)/filtered,
           lost_merged = (denoisedR-merged)/denoisedR, # Proportion of reads lost to merging
           prop_chimera = (merged-nonchim)/merged) %>% 
    pivot_longer(where(is.numeric), names_to = 'variable', values_to = 'values') 
}


###############
### PLOTS ######
#################

plot_track_change <- function(track_change) {
  change_vars <- c('prop_chimera', 'lost_merged', 'lost_noise', 'lost_filt', 'lost_Ns')
  track_change %>% 
    filter(variable %in% change_vars) %>% 
    mutate(variable = factor(variable, level = change_vars)) %>% 
    ggplot(aes(y = variable, x = values)) +
    geom_jitter() + theme_minimal() # overall very low chimeric rate
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

