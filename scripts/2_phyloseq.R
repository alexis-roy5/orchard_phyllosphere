library(pacman)
p_load(tidyverse, phyloseq, magrittr, decontam, Biostrings, readxl)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/phyloseq_functions.R")


##########
# Setup ###
############

setwd("/Users/alexisroy/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere/")
path.out.ps <- "./2023/out/"

# Metadata
meta_raw <- read_excel("./2023/data/ITS/4_taxonomy_ITS_real/Verger_2023_1.xlsx", sheet = "Metadata") # change ME to your Metadata sheet/file
meta <- meta_raw %>%
  mutate(unique = str_remove(unique, "-16S" ),
         TreeID = str_extract(unique, "(...)-(..)-")) %>% 
  unite("TreeID", TreeID, replicate, sep = "", remove = FALSE)


# adding DNA concetration to the meatadata
DNA_meta_raw <- read_excel("./2023/data/ITS/4_taxonomy_ITS_real/Verger_2023_1.xlsx", sheet = "Biomol")

DNA_meta <- DNA_meta_raw %>%
  select(unique, Concentration = `ng/ul`) %>%
  as.data.frame() 

samples_meta <- left_join(meta, DNA_meta, by = "unique")


#########
# 16S ####
###########

path_16S <- './2023/data/16S/4_taxonomy_16S'
taxa_16S_genus <- read_rds(file.path(path_16S, 'taxonomy.rds'))
taxa_16S_species <- read_rds(file.path(path_16S, 'taxonomy_species.rds'))
seqtab_16S <- read_rds(file.path(path_16S, 'seqtab.rds'))

# Somehow assignSpecies enlève tous les autres rangs autres que Genre et Espèce. 
# Fixons cela!
Species_16S <- taxa_16S_species[,2]
names(Species_16S) <- rownames(taxa_16S_species)
taxa_16S <- cbind(taxa_16S_genus, Species_16S) %>% data.frame

# Keep samples with metadata info
seqtab_16S_sam <- subset_samples(seqtab_16S, sample.names.phylo)

# maximum de hits de ASVs dans le control négatif
max_value <- max(seqtab_16S_ctrl[1, ])
print(max_value)

# Subset ASVs
taxa_16S_sam <- subset_asvs(taxa_16S, seqtab_16S_sam, 100) 

# looking at the repartition of sample size - to choose the treshold for "remove_ultra_rare"
rowSums(seqtab_16S_sam) %>% 
  sort() %>% 
  hist(breaks = 50, xlab = "sample size", ylab = "nb of samples", xaxt = "n", main = "")  

axis(1, at = pretty(rowSums(seqtab_16S_sam), n = 20))  # adding ticks 


# Remove near-empty samples - LOOK at the number
seqtab_16S_sam_filt <- remove_ultra_rare(seqtab_16S_sam, taxa_16S_sam, 4000)
dim(seqtab_16S_sam); dim(seqtab_16S_sam_filt); dim(taxa_16S_sam)


# Finding the sample that are near-empty 
setdiff(rownames(seqtab_16S_sam), rownames(seqtab_16S_sam_filt)) 

# Phyloseq object
ps_16S <- phyloseq(
  tax_table(taxa_16S_sam),
  otu_table(seqtab_16S_sam_filt, taxa_are_rows = FALSE),
  sample_data(samples)
)

saveRDS(ps_16S, paste0(path.out.ps, "ps_16S.rds"))

#########
# ITS ####
###########

path_ITS <- './2023/data/ITS/'
taxa_ITS <- read_rds(file.path(path_ITS, '4_taxonomy_ITS_real/taxonomy.RDS'))
seqtab_ITS <- read_rds(file.path(path_ITS, '4_taxonomy_ITS_real/seqtab.RDS'))

# add sequencing depth to sample data
samples_meta_updated <- seqtab_ITS %>% 
  rowSums() %>% 
  data.frame(seqDepth = .) %>% 
  rownames_to_column(var = "sample") %>% 
  left_join(samples_meta,. ,by = "sample") %>% column_to_rownames(var = "sample")

# Keep samples with metadata info
seqtab_ITS_sam <- subset_samples(seqtab_ITS, samples_meta$sample) # s'assurer que les infos de noms sont identiques pour 16S, ITS, trnl, etc.

# Subset ASVs (ex: 100 hits)
taxa_ITS_sam <- subset_asvs(taxa_ITS, seqtab_ITS_sam, 10)

# looking at the repartition of sample size - to choose the treshold for "remove_ultra_rare"
rowSums(seqtab_ITS_sam) %>% 
  sort() %>% 
  hist(breaks = 100, xlab = "sample size", ylab = "nb of samples", xaxt = "n", main = "")  

axis(1, at = pretty(rowSums(seqtab_ITS_sam), n = 20))  # adding ticks 

# Remove near-empty samples - LOOK at the number and change it in consequences
seqtab_ITS_sam_filt <- remove_ultra_rare(seqtab_ITS_sam, taxa_ITS_sam, 2000) 

dim(seqtab_ITS_sam); dim(seqtab_ITS_sam_filt); dim(taxa_ITS_sam)

# Finding the sample that are near-empty 
setdiff(rownames(seqtab_ITS_sam), rownames(seqtab_ITS_sam_filt))

# Phyloseq objects
ps_ITS <- phyloseq(
  tax_table(taxa_ITS_sam),
  otu_table(seqtab_ITS_sam_filt, taxa_are_rows = FALSE),
  sample_data(samples_meta_updated)
)

saveRDS(ps_ITS, paste0(path.out.ps, "ps_ITS.rds"))
