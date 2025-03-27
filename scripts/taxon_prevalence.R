# Compute taxon prevalence, top 5 genera
library(pacman)
p_load(tidyverse, magrittr)

(ps <- readRDS('2023/out/ps_16S.rds'))
start.time <- Sys.time()
melted <- psmelt(ps)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
melted_fast <- psmelt_dt(ps)
end.time <- Sys.time()
end.time - start.time


# Most prevalent genera
melted %>% 
  select(Sample, Abundance, Genus) %>% 
  filter(Abundance>0) %>% 
  group_by(Sample, Genus) %>% 
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>% 
  group_by(Genus) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))

# vs. Genera of most prevalent ASVs ? 
  
  
