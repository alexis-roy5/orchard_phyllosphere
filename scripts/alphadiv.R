library(pacman)
p_load(tidyverse, phyloseq)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')

ps_ITS <- readRDS("2023/out/ps_ITS.rds")
ps = rarefy_even_depth(ps_ITS,  
                       sample.size=0.95*min(sample_sums(ps_ITS)), 
                       replace=F, rngseed = 42)

samdat <- ps %>% 
  samdat_as_tibble %>% 
  select(Sample, time,site, cultivar, type, replicate, practice, susceptibility, date) %>% 
  mutate(treeID = paste0(site,'_',cultivar,'_', replicate))

map(c('Shannon', 'Simpson', 'Tail'), function(idx){
  estimate_diversity(ps, index = idx) %>% 
    as.data.frame() %>% 
    setNames(idx)
}) %>% list_cbind %>% 
  rownames_to_column('Sample') %>% 
  tibble %>% 
  left_join(samdat, by = 'Sample') %>% 
  write_tsv('2023/data/orchard_ITS_div_2023.tsv')
