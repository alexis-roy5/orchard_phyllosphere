## phyloseq pipeline ##
library(pacman)
p_load(tidyverse, phyloseq, magrittr, kableExtra, Biostrings, readxl)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/phyloseq_functions.R")


##########
# Setup ###
############

setwd("/Users/alexisroy/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere/")
path.out.ps <- "./2023/out/" # place for the final object phyloseq

# metadata
meta_raw <- read_excel("./2023/data/ITS/4_taxonomy_ITS/Verger_2023_1.xlsx", sheet = "Metadata") # change ME to your Metadata sheet/file
meta <- meta_raw %>%
  mutate(unique = str_remove(unique, "-16S" ),
         TreeID = str_extract(unique, "(...)-(..)-")) %>% 
  unite("TreeID", TreeID, replicate, sep = "", remove = FALSE)

#########
# ITS ####
###########

path_ITS <- './2023/data/ITS/'
taxa_ITS <- read_rds(file.path(path_ITS, '4_taxonomy_ITS/taxonomy.RDS'))
seqtab_ITS <- read_rds(file.path(path_ITS, '4_taxonomy_ITS/seqtab.RDS'))

# add sequencing depth to sample data
samples_meta_updated <- seqtab_ITS %>% 
  rowSums() %>% 
  data.frame(seqDepth = .) %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(meta,. ,by = "sample") %>% column_to_rownames(var = "sample")

# keep samples with metadata info
seqtab_ITS_sam <- subset_samples(seqtab_ITS, meta$sample) # s'assurer que les infos de noms sont identiques pour 16S, ITS, trnl, etc.

# subset ASVs (ex: 10 hits)
taxa_ITS_sam <- subset_asvs(taxa_ITS, seqtab_ITS_sam, 10)

# looking at the repartition of sample size - to choose the treshold for "remove_ultra_rare"
rowSums(seqtab_ITS_sam) %>% 
  sort() %>% 
  hist(breaks = 100, xlab = "sample size", ylab = "nb of samples", xaxt = "n", main = "")  

axis(1, at = pretty(rowSums(seqtab_ITS_sam), n = 20))  # adding ticks 

# remove near-empty samples - LOOK at the number and change it in consequences
seqtab_ITS_sam_filt <- remove_ultra_rare(seqtab_ITS_sam, taxa_ITS_sam, 2000) 

dim(seqtab_ITS_sam); dim(seqtab_ITS_sam_filt); dim(taxa_ITS_sam)

# finding the sample that are near-empty 
setdiff(rownames(seqtab_ITS_sam), rownames(seqtab_ITS_sam_filt))

# phyloseq objects
ps_ITS <- phyloseq(
  tax_table(taxa_ITS_sam),
  otu_table(seqtab_ITS_sam_filt, taxa_are_rows = FALSE),
  sample_data(samples_meta_updated)
)

saveRDS(ps_ITS, paste0(path.out.ps, "ps_ITS.rds"))


#################################
### STATS for phyloseq object ###
#################################
# do the same thing after for the non rarefied object 
ps.rarefied.ITS = rarefy_even_depth(ps_ITS, rngseed=1, sample.size=2500, replace=F) # rarefying at 2500
saveRDS(ps.rarefied.ITS, paste0(path.out.ps, "ps_ITS_rarefied.rds"))

################
# ASV.SEQUENCE #
################
# goal: prepare a table with information of ASV and sequences
# prep data
ps.stats.prep <- function(ps, barcode) {
  asv <- ps %>% otu_table 
  seq_per_sam <- rowSums(asv)
  asv_per_sam <- rowSums(asv>0)
  asv_prevalence <- colSums(asv>0) # nb ASV se retrouve dans au moins N sample
  num_sam <- nrow(asv)
  tibble(
    Dataset = barcode,
    Seq = sum(asv),
    ASVs = ncol(asv),
    N = num_sam,
    Mean_seq = mean(seq_per_sam),
    SD_seq = sd(seq_per_sam),
    Min_seq = min(seq_per_sam),
    Max_seq = max(seq_per_sam),
    Mean_asv = mean(asv_per_sam),
    SD_asv = sd(asv_per_sam),
    Min_asv = min(asv_per_sam),
    Max_asv = max(asv_per_sam),
    Mean_prev = mean(asv_prevalence),
    SD_prev = sd(asv_prevalence),
    Min_prev = min(asv_prevalence),
    Max_prev = max(asv_prevalence)
  )
} 

# execute computations
ps.stats <- ps.stats.prep(ps.rarefied.ITS, "ITS") %>% 
  mutate(across(where(is.numeric), ~ format(round(., 0),big.mark=',')))

ps.stats.k <- kable(ps.stats, "html", align = "c") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c( # format table
    "Dataset" = 1,
    "Sequences" = 1,
    "ASVs" = 1,
    "Samples" = 1,
    "Mean ± SD" = 2, 
    "[Min, Max]" = 2, 
    "Mean ± SD" = 2, 
    "[Min, Max]" = 2, 
    "Mean ± SD" = 2, 
    "[Min, Max]" = 2
  )) %>%  
  add_header_above(c(
    " " = 4, 
    "Sequences per sample" = 4, 
    "ASVs per sample" = 4, 
    "ASV prevalence" = 4
  )) %>%
  row_spec(0, extra_css = "display: none;")  # Hide the original column names


html_file <- file.path("./2023/out/table_asv_sequence_rarefaction.html") # path for table
save_kable(ps.stats.k, file = html_file) # save
