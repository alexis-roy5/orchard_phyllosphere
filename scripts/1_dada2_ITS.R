#ITS
#  ml StdEnv/2023 r/4.4.0 mugqic/cutadapt/2.10
setwd("/home/def-ilafores/analysis/orchard_phyllosphere")

library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel)
source('./scripts/myFunctions.R')
source('./scripts/mergePairsRescue.R')



# CONFIG
barcode <- '-ITS'
prefix <- '2023-'
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  # Is it the good ones ?
REV <- "GCTGCGTTCTTCATCGATGC"

ncores <- 48
path_data <- paste0('./2023/data/',barcode) # change the YEAR
path_raw <- paste0(path_data, '/0_raw')
if(!dir.exists(path_raw)) message("create.directory.'/0_raw'")

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE)) # change ME to your path for RAW data
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(fnFs, get.sample.name, prefix = year, USE.NAMES = FALSE)
# write_delim(data.frame(sample.names), paste0('data/sample_names_',barcode,'.tsv'))

########################
# 1. N-FILTERING ########
##########################

### PRE-FILTER (no qc, just remove Ns)
fnFs.filtN <- file.path(path_data, "1_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path_data, "1_filtN", basename(fnRs))
out.N <- filterAndTrim(fnFs, fnFs.filtN,
                         fnRs, fnRs.filtN, 
                         rm.lowcomplex = TRUE, # added because of https://github.com/benjjneb/dada2/issues/2045#issuecomment-2452299127
                         maxN = 0, 
                         multithread = ncores)


head(out.N)

###########################
# 2. PRIMER REMOVAL ########
#############################

# Analyse primer occurence
primer_occurence(fnFs.filtN, fnRs.filtN, FWD, REV)

### CUTADAPT
#cutadapt <-"/Users/alexisroy/miniconda3/envs/cutadapt/bin/cutadapt" # my computer
cutadapt <- '/cvmfs/soft.mugqic/CentOS6/software/cutadapt/cutadapt-2.10/bin/cutadapt' # IP34. CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path_data, "2_cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt multicore (each sample is run single-core)
mclapply(seq_along(fnFs), run_cutadapt, mc.cores = ncores)

# Check if it worked?
primer_occurence(fnFs.cut, fnRs.cut, FWD, REV)

##############################
# 3. QUALITY FILTERING ########
################################

### Filter and trim for real
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Check quality
plotQualityProfile(cutFs[10:21])
plotQualityProfile(cutRs[10:21])

# Filter samples; define out files
filtFs <- file.path(path_data, "3_filtered", basename(cutFs))
filtRs <- file.path(path_data, "3_filtered", basename(cutRs))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxEE = c(2, 2), 
                     truncQ = 2,
                     minLen = 100,
                     rm.phix = TRUE, 
                     compress = TRUE, 
                     multithread = ncores) 
plotQualityProfile(filtFs[10:21])
plotQualityProfile(filtRs[10:21])

####################################
# 4. ERROR MODEL AND CLEANUP ########
######################################

# Filtering with minlen may yield empty samples (e.g. neg. controls);
# list files that did survive filtering:
filtFs_survived <- filtFs[file.exists(filtFs)]
filtRs_survived <- filtRs[file.exists(filtRs)]


# list files that did NOT survived
names(filtFs[!file.exists(filtFs)])
names(filtRs[!file.exists(filtRs)])

# Learn errors from the data
errF <- learnErrors(filtFs_survived, multithread = ncores)
errR <- learnErrors(filtRs_survived, multithread = ncores)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Infer sample composition
dadaFs <- dada(filtFs_survived, err = errF, 
               pool = 'pseudo', multithread = ncores)
dadaRs <- dada(filtRs_survived, err = errR, 
               pool = 'pseudo', multithread = ncores) # pseudo pool. putting all the error "rates" information together but a heuristic.


####################################
######## MERGING TEST ##############
####################################

# only using the forward reads
forwardFs<- dadaFs

# mergePairs()
merged <- mergePairs(dadaFs, filtFs_survived, dadaRs, filtRs_survived, verbose=TRUE)

# mergePairs(..., just_Concatenate)
concatenated <- mergePairs(dadaFs, filtFs_survived, dadaRs, filtRs_survived, justConcatenate = TRUE, verbose=TRUE)

# Modified version of mergePairs that rescues non-merged reads by concatenation
# source('scripts/mergePairsRescue.R')
rescued <- mergePairsRescue(
  dadaFs, filtFs_survived, 
  dadaRs, filtRs_survived,
  returnRejects = TRUE,
  minOverlap = 12,
  maxMismatch = 0,
  rescueUnmerged = TRUE
)

# Intersect the merge and concat; allows merge to fail when overlap is mismatched,
# but recovers non-overlapping pairs by concatenating them. 
# Motivated by https://github.com/benjjneb/dada2/issues/537#issuecomment-412530338

path.tax <- file.path(path_data, "4_taxonomy_ITS")
if(!dir.exists(path.tax)) dir.create(path.tax)

# seqtab building
# Remove chimeras
# WRITE OUT without chimeras

### ASSIGN TAXONOMY
# avoir sa database dans le dossier "reference_database"
reference_path <- paste0(path_data, '/reference_database/sh_general_release_dynamic_s_04.04.2024.fasta') # change to your reference database

seqtab.nochim.Fs <- write.seqtab.taxonomic(forwardFs,path.tax, reference_path, ncores)
seqtab.nochim.merged <- write.seqtab.taxonomic(merged, path.tax, reference_path, ncores)
seqtab.nochim.concatenated <- write.seqtab.taxonomic(concatenated, path.tax, reference_path, ncores)
seqtab.nochim.rescued <- write.seqtab.taxonomic(rescued, path.tax, reference_path, ncores)

### TRACK PIPELINE READS 

# I don't find the need for track_change in this experiment ("merge test")

track_change <- track_dada(out.N = out.N, out = out,
                           dadaFs = dadaFs, dadaRs = dadaRs,
                           mergers = mergers_pooled,
                           seqtab.nochim = seqtab.nochim)

path.out <- file.path("./2023", "out")
if(!dir.exists(path.out)) dir.create(path.out)

track_change %>%
  plot_track_change() %>% 
  ggsave(paste0('./2023/out/change_',barcode,'.pdf'), plot = ., 
         bg = 'white', width = 1600, height = 1200, 
         units = 'px', dpi = 180)

################################
#### TAXONOMIC ASSIGNATION ####
###############################

## optimiser le début de mon script et le faire pour les 4 méthodes

tax_forward <- read_rds(file.path(path.tax, "taxonomy.forwardFs.RDS"))
tax_merged <- read_rds(file.path(path.tax, "taxonomy.merged.RDS"))
tax_concatenated <- read_rds(file.path(path.tax, "taxonomy.concatenated.RDS"))
tax_rescued <- read_rds(file.path(path.tax, "taxonomy.rescued.RDS"))

seq_forward <- read_rds(file.path(path.tax, "seqtab.forwardFs.RDS"))
seq_merged <- read_rds(file.path(path.tax, "seqtab.merged.RDS"))
seq_concatenated <- read_rds(file.path(path.tax, "seqtab.concatenated.RDS"))
seq_rescued <- read_rds(file.path(path.tax, "seqtab.rescued.RDS"))


# SETUP

tax_tables <- list(
  forward = tax_forward,
  merged = tax_merged,
  concatenated = tax_concatenated,
  rescued = tax_rescued
)

seq_tables <- list(
  forward = seq_forward,
  merged = seq_merged,
  concatenated = seq_concatenated,
  rescued = seq_rescued
)
# function: transform a dataframe into a longer one
transform_df_long <- function(object) { 
  as.data.frame(object) %>% # transform an object/vector into a data.frame 
    rownames_to_column('Sample') %>%
    pivot_longer(., 
                 cols = -Sample,  # exclude the 'Sample' column
                 names_to = "ASV", 
                 values_to = "Abundance") %>% 
  filter(Abundance > 0)
}

transform_df_long(seq_forward) %>% View()

setup <- function(taxonomic_table, seqtab) {
  seq_long <- transform_df_long(seqtab) # changing seqtab
  tax_modif <- as.data.frame(rownames_to_column(taxonomic_table, var = "ASV")) # changing tax table
  seq_RA <- seq_long %>% 
    group_by(Sample) %>% 
    mutate(Relative_Abundance = Abundance / sum(Abundance)) %>% 
    left_join(tax_modif, by = "ASV")
  return(seq_RA) 
}


results <- map2(tax_tables, seq_tables, setup)  # contains all the dataframe in 1 list
# map ; returns a list of outputs ; with 2 inputs

# list of assignations & options 
assignation_levels_imp <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
options_methods <- c("forward", "merged", "concatenated", "rescued")

# function to know the pourcentage of assignation of one level
pourcentage_assignation.1 <- function(taxonomic_table, level, weighted_by_RA = FALSE) {
  result <- taxonomic_table %>%
    group_by(Sample) %>%  
    
    { if (!weighted_by_RA) # calculate percentage without weighting by relative abundance
      summarise(., Pourcentage = sum(!is.na(.data[[level]])) / n() * 100, .groups = "drop") #.data[[level]]) : input of pipe = .data, $level dynamic = [[level]]
      else summarise(., Pourcentage = sum(Relative_Abundance * !is.na(.data[[level]])) * 100, .groups = "drop")} %>% # Calculate percentage weighted by relative abundance
    
    mutate(Level = level)  # keep track of level information
  return(result)
}

pourcentage_assignation(results$forward, "Genus")
pourcentage_assignation.1(results$forward, "Genus")

# function to know the pourcentage of assignation of ALL levels for one table
all_assignation <- function(taxonomic_table, levels_list, weighted_by_RA = FALSE) {
  map_dfr(.x = levels_list, ~ pourcentage_assignation(taxonomic_table, .x, weighted_by_RA)) # applying pourcentage assignation for each table
}

# for ALL the taxonomic tables
frame_assignation <- function(results, assignation_levels, weighted_by_RA = FALSE) {
  map_dfr(names(results), function(table_name) { # function anonyme, like a FOR loop
    taxonomic_table_real <- results[[table_name]]  # for each table_name
    all_assignation(taxonomic_table = taxonomic_table_real, level = assignation_levels, weighted_by_RA) %>%
      mutate(TaxTable = table_name)  # Add table name as a column
  })
}

df.sequence <- frame_assignation(results, assignation_levels_imp, weighted_by_RA = FALSE)
df.weighted.RA <- frame_assignation(results, assignation_levels_imp, weighted_by_RA = TRUE)
# map_dfr : Binds results by row (bind_rows()) ; with 1 input

##########
## PLOT ##
##########

library(ggridges)

# function: create a plot on the level
create_rich_plot_unique <- function(data, level) {
  # filter the data based on the provided level
  filtered_data <- data %>% 
    filter(Level == level) %>% 
    mutate(TaxTable = factor(TaxTable, 
                             levels = c("forward", "merged", "concatenated", "rescued")))
  
  # creating the plot
  ggplot(filtered_data, aes(x = Pourcentage, y = TaxTable, fill = TaxTable)) +
    labs(title = level) +
    geom_density_ridges2(
      alpha = 0.6, show.legend = FALSE, 
      scale = 3 #, stat = "binline", bins = 100  # histogram-like ridges
    ) +
    labs(x = 'Sequence assignation percentage', y = 'Methods') +
    theme_classic() +
    scale_x_continuous(limits = c(0, 100))
}

# Assuming df.weighted.RA is your data frame

create_rich_plot <- function(data_frame) {
ggplot(data_frame %>% 
         mutate(Level =factor(Level, levels = rev(c("Phylum", "Class", "Order", "Family", "Genus", "Species")))),
                aes(x = Pourcentage, y = Level, fill = TaxTable)) +
  geom_density_ridges(scale = 1, alpha = 0.5) +
  theme_ridges() +
  labs(x = "Percentage", y = "Tax Rank", fill = "Merging options") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(50, 100))
}

sequence.plot <- create_rich_plot(df.sequence)
sequence.weighted.RA.plot <- create_rich_plot(df.weighted.RA)

###################
### MEAN and SD ###
###################
add.mean.sd <- function(data.frame) {
  # specified levels
  data.frame <- data.frame %>%
    mutate(Level = factor(Level, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")))
  
  # group by TaxTable and Level, then calculate mean and SD
  mean.sd <- data.frame %>% 
    group_by(TaxTable, Level) %>% 
    summarise(
      Mean = mean(Pourcentage),  
      Standard_deviation = sd(Pourcentage),
      .groups = "drop"
    )
  
  return(mean.sd)

}

Sequence.mean.sd <- add.mean.sd(df.sequence)
ASV.mean.sd <- add.mean.sd(df.ASV)

