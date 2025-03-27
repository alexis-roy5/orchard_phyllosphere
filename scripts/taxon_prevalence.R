# Compute taxon prevalence, top 5 genera
library(pacman)
p_load(tidyverse, magrittr)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")

(ps <- readRDS('2023/out/ps_ITS.rds'))

melted  <- psflashmelt(ps)

# Most prevalent Genera
melted %>%
  filter(practice == 'Organic') %>% 
  select(Sample, Abundance, Genus) %>% 
    group_by(Sample, Genus) %>% 
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>% 
    group_by(Genus) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))

# vs. Genera of most prevalent ASVs ? 
asv_prevalence_table <- function(filtered_table, taxRank, num = Inf){
  filtered_table %>% 
    select(Sample, Abundance, !!sym(taxRank), OTU) %>% 
    group_by(OTU, !!sym(taxRank)) %>% 
    summarise(n = n(), .groups = 'drop') %>%
    arrange(desc(n)) %>% 
    select(Genus, n, OTU) %>% 
    head(n = num) %>% 
    mutate(ASV_rank = row_number())
}

asv_org <- melted %>% 
  filter(practice == 'Organic') %>% 
  asv_prevalence_table('Genus')
  
asv_conv <- melted %>% 
  filter(practice == 'Conventional') %>% 
  asv_prevalence_table('Genus') 

joined_table <- full_join(
  asv_conv, asv_org, 
  by = c('OTU', 'Genus'),
  suffix = c('.conv', '.org')) %>%
  rename(ASV = OTU) %>% 
  mutate(ASV_label = paste0("ASV_", row_number()))

  
joined_table %>% 
  select(-ASV) %>% 
  filter(ASV_rank.org %in% seq(1:10) |
           ASV_rank.conv %in% seq(1:10)) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling()

joined_table %>% 
  select(ASV_label, ASV) %>% 
  write_tsv('2023/out/ASV_labels.tsv')
