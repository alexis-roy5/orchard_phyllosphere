#ITS
#  ml StdEnv/2023 r/4.4.0 mugqic/cutadapt/2.10
library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel)
source('scripts/myFunctions.R')
source('scripts/mergePairsRescue.R')

setwd("/home/def-ilafores/analysis/orchard_phyllosphere")

# CONFIG
barcode <- 'ITS'
FWD <- "CTTGGTCATTTAGAGGAAGTAA" 
REV <- "GCTGCGTTCTTCATCGATGC"

ncores <- 24
path_data <- paste0('./2023/data/ITS',barcode)
path_raw <- paste0(path_data, '/0_raw')
if(!dir.exists("0_raw")) message("create.directory.'/0_raw'")

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

(sample.names <- sapply(fnFs, get.sample.name, USE.NAMES = FALSE))
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
#head(out.N)

###########################
# 2. PRIMER REMOVAL ########
#############################

# Analyse primer occurence
primer_occurence(fnFs.filtN, fnRs.filtN, FWD, REV)

### CUTADAPT
#cutadapt <- "/Users/jorondo/miniconda3/envs/cutadapt/bin/cutadapt" 
# cutadapt <- '/cvmfs/soft.mugqic/CentOS6/software/cutadapt/cutadapt-2.10/bin/cutadapt' # CHANGE ME to the cutadapt path on your machine
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

# Learn errors from the data
errF <- learnErrors(filtFs_survived, multithread = ncores)
errR <- learnErrors(filtRs_survived, multithread = ncores)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Infer sample composition
dadaFs <- dada(filtFs_survived, err = errF, 
               pool = 'pseudo', multithread = ncores)
dadaRs <- dada(filtRs_survived, err = errR, 
               pool = 'pseudo', multithread = ncores)

# Modified version of mergePairs that rescues non-merged reads by concatenation
# source('scripts/mergePairsRescue.R')
mergers_pooled <- mergePairsRescue(
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
path.tax <- file.path(path_data, "4_taxonomy")
if(!dir.exists(path.tax)) dir.create(path.tax)

seqtab <- makeSequenceTable(mergers_pooled) # makeSequenceTable(dadaFs) ## to use FWD READS ONLY

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread = ncores, verbose = TRUE
)

dim(seqtab); dim(seqtab.nochim) # from 11139 to 9043
# WRITE OUT
write_rds(seqtab.nochim, paste0(path.tax,'/seqtab.RDS'))

### TRACK PIPELINE READS
track_change <- track_dada(out.N = out.N, out = out,
                           dadaFs = dadaFs, dadaRs = dadaRs,
                           mergers = mergers_pooled,
                           seqtab.nochim = seqtab.nochim)

track_change %>% 
  filter(values>=0) %>% 
  plot_track_change() %>% 
  ggsave(paste0('out/change_ITS',barcode,'.pdf'), plot = ., 
         bg = 'white', width = 1600, height = 1200, 
         units = 'px', dpi = 180)

### ASSIGN TAXONOMY
taxa <- assignTaxonomy(
  seqtab.nochim, 
  paste0(path_data, '/sh_general_release_dynamic_04.04.2024_dev.fasta'), 
  multithread= min(ncores, 72), tryRC = TRUE, verbose = TRUE)

taxa_fixed <- taxa %>% as.data.frame() %>%
  mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .))) 

write_rds(taxa_fixe, paste0(path.tax,'/taxonomy.RDS'))

