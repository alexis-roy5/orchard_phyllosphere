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


meta_samples <- meta %>% meta_samples <- meta %>% replicate
  filter(time != "None" & control == 'FALSE')

meta_ctrl <- meta %>% 
  filter(time =="None" | control=='TRUE')
  
sample.names.phylo <- meta_samples$unique
ctrl.names <- meta_ctrl$unique



# adding DNA concetration to the meatadata
  DNA_meta_raw <- read_excel("./2023/data/ITS/4_taxonomy_ITS_real/Verger_2023_1.xlsx", sheet = "Biomol")
  
  DNA_meta <- DNA_meta_raw %>%
    select(unique, Concentration = `ng/ul`) %>%
    as.data.frame() 
  
samples <- left_join(meta, DNA_meta, by = "unique") %>% 
  column_to_rownames("unique") 
  
  
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
seqtab_16S_ctrl <- subset_samples(seqtab_16S, ctrl.names)

# maximum de hits de ASVs dans le control négatif
max_value <- max(seqtab_16S_ctrl[1, ])
print(max_value)

# Subset ASVs
taxa_16S_sam <- subset_asvs(taxa_16S, seqtab_16S_sam, 100) 
taxa_16S_ctrl <- subset_asvs(taxa_16S, seqtab_16S_ctrl, 100) 

# looking at the repartition of sample size - to choose the treshold for "remove_ultra_rare"
rowSums(seqtab_16S_sam) %>% 
  sort() %>% 
  hist(breaks = 50, xlab = "sample size", ylab = "nb of samples", xaxt = "n", main = "")  

axis(1, at = pretty(rowSums(seqtab_16S_sam), n = 20))  # adding ticks 


# Remove near-empty samples - LOOK at the number
seqtab_16S_sam_filt <- remove_ultra_rare(seqtab_16S_sam, taxa_16S_sam, 4000)
seqtab_16S_ctrl_filt <- remove_ultra_rare(seqtab_16S_ctrl, taxa_16S_ctrl, 4000) 
dim(seqtab_16S_sam); dim(seqtab_16S_sam_filt); dim(taxa_16S_sam)
dim(seqtab_16S_ctrl); dim(seqtab_16S_ctrl_filt); dim(taxa_16S_ctrl)

  # Finding the sample that are near-empty 
    setdiff(rownames(seqtab_16S_sam), rownames(seqtab_16S_sam_filt)) 

# Phyloseq object
ps_16S <- phyloseq(
  tax_table(taxa_16S_sam),
  otu_table(seqtab_16S_sam_filt, taxa_are_rows = FALSE),
  sample_data(samples)
)

saveRDS(ps_16S, paste0(path.out.ps, "ps_16S.rds"))

##########################################################
##########################################################

ps_16S_ctrl <- phyloseq(
  tax_table(taxa_16S_ctrl),
  
  otu_table(seqtab_16S_ctrl_filt, taxa_are_rows = FALSE),
  sample_data(samples)
)

# Export asvs as fasta
asv_to_fasta(seqtab_16S_sam_filt, file.path(path_16S, '4_taxonomy/asv.fa'))

#########
# ITS ####
###########

path_ITS <- './2023/data/ITS/'
taxa_ITS <- read_rds(file.path(path_ITS, '4_taxonomy_ITS_real/taxonomy.RDS'))
seqtab_ITS <- read_rds(file.path(path_ITS, '4_taxonomy_ITS_real/seqtab.RDS'))

rownames(seqtab_ITS) <- paste0("2023-", rownames(seqtab_ITS)) # adding the prefix '2023_' because it was missing from the start (raw_data) )

# Keep samples with metadata info
seqtab_ITS_sam <- subset_samples(seqtab_ITS, sample.names.phylo) # s'assurer que les infos de noms sont identiques pour 16S, ITS, trnl, etc.
seqtab_ITS_ctrl <- subset_samples(seqtab_ITS, ctrl.names)

# Subset ASVs (ex: 100 hits)
taxa_ITS_sam <- subset_asvs(taxa_ITS, seqtab_ITS_sam, 10)
taxa_ITS_ctrl <- subset_asvs(taxa_ITS, seqtab_ITS_ctrl, 100) # ne fonctionne pas avec 1 sample 

# looking at the repartition of sample size - to choose the treshold for "remove_ultra_rare"
rowSums(seqtab_ITS_sam) %>% 
  sort() %>% 
  hist(breaks = 100, xlab = "sample size", ylab = "nb of samples", xaxt = "n", main = "")  

axis(1, at = pretty(rowSums(seqtab_ITS_sam), n = 20))  # adding ticks 

# Remove near-empty samples - LOOK at the number and change it in consequences
seqtab_ITS_sam_filt <- remove_ultra_rare(seqtab_ITS_sam, taxa_ITS_sam, 2000) 
seqtab_ITS_ctrl_filt <- remove_ultra_rare(seqtab_ITS_ctrl, taxa_ITS_ctrl, 2000)

dim(seqtab_ITS_sam); dim(seqtab_ITS_sam_filt); dim(taxa_ITS_sam)
dim(seqtab_ITS_ctrl); dim(seqtab_ITS_ctrl_filt); dim(taxa_ITS_ctrl)

# Finding the sample that are near-empty 
setdiff(rownames(seqtab_ITS_sam), rownames(seqtab_ITS_sam_filt))

# Phyloseq objects
ps_ITS <- phyloseq(
  tax_table(taxa_ITS_sam),
  otu_table(seqtab_ITS_sam_filt, taxa_are_rows = FALSE),
  sample_data(samples)
  )

saveRDS(ps_ITS, paste0(path.out.ps, "ps_ITS.rds"))

ps_ITS_ctrl <- phyloseq(
  tax_table(taxa_ITS_ctrl),
  otu_table(seqtab_ITS_ctrl_filt, taxa_are_rows = FALSE),
  sample_data(meta_ctrl_ITS)
)

# Export asvs as fasta
asv_to_fasta(seqtab_ITS_sam_filt, file.path(path_ITS, '4_taxonomy/asv.fa'))

# //DEV
# Decontamination DECONTAM
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")

# code Alex
# prune rows NA in my phyloseq object - make it numeric
ps_16S_clean <- prune_samples(!is.na(sample_data(ps_16S)$Concentration), ps_16S)
sample_data(ps_16S_clean)$Concentration <- as.numeric(sample_data(ps_16S_clean)$Concentration)
## j'ai enelever des échantillons, est-ce correct ?

contamdf.freq <- isContaminant(ps_16S_clean, method="frequency", conc="Concentration") # conc needs to be grater than 0
View(contamdf.freq) # $contaminant=TRUE if $p < 0.1

# how many contaminants
table(contamdf.freq$contaminant) # TRUE = contaminants 

list(which(contamdf.freq$contaminant)) # la table est en fréquence d'ASV décroissant, montre leur importance
                                        # exemple: 3 = l'ASV le 3e plus fréquent

plot_frequency(ps_16S_clean, taxa_names(ps_16S_clean)[c(1,119)], conc="Concentration") + 
  xlab("DNA Concentration (ng/ul)")
# ça ressemble pas à ce qu'il dit 

#############################
# code de Jo

ps_16S_load <- ps_16S %>%
  prune_samples(samples(.)$Concentration>0, .)

contam_freq<- ps_16S_load %>% 
  isContaminant(method = 'frequency', conc = 'Concentration')

table(contam_freq$contaminant)
head(which(contam_freq$contaminant))

ps_16S_load %>% 
  plot_frequency(., taxa_names(.)[c(516,1128)], conc="bacterial_load") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

