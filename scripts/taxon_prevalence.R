# Compute taxon prevalence, top 5 genera
library(pacman)
p_load(tidyverse, magrittr, kableExtra)
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
asv_prevalence_table <- function(filtered_table, taxRank, num = Inf) {
  filtered_table %>% 
    select(Sample, Abundance, !!sym(taxRank), OTU) %>% 
    group_by(OTU, !!sym(taxRank)) %>% 
    summarise(Prev = n(), .groups = 'drop') %>%
    arrange(desc(Prev)) %>% 
    select(!!sym(taxRank), Prev, OTU) %>% 
    head(n = num) %>% 
    mutate(ASV_rank = row_number())
}

melted_org <- melted %>% 
  filter(practice == 'Organic')

melted_conv <- melted %>% 
  filter(practice == 'Conventional')

asv_org <- asv_prevalence_table(melted_org, 'Genus')
asv_conv <- asv_prevalence_table(melted_conv, 'Genus')

joined_table <- full_join(
  asv_conv, asv_org, 
  by = c('OTU', 'Genus'),
  suffix = c('.conv', '.org')) %>%
  rename(ASV = OTU) %>% 
  mutate(ASV_label = paste0("ASV_", row_number()))

# Comparative table  
## Header names set dynamically to include sample count:
nsam_org <- melted_org %>% pull(Sample) %>% unique %>% length
nsam_conv <- melted_conv %>% pull(Sample) %>% unique %>% length

top_headers_names <- setNames(
  c(1, 2, 2, 1), 
  c(" ", 
    paste0("Conventional (n=", nsam_conv, ")"), 
    paste0("Organic (n=", nsam_org, ")" ),
    " ")
)

# Table :
joined_table %>% 
  mutate(compound_score = ASV_rank.conv + ASV_rank.org,
         Prev.org = paste0(Prev.org, " (",round(100*Prev.org/nsam_org,0), " %)"),
         Prev.conv = paste0(Prev.conv, " (",round(100*Prev.conv/nsam_conv,0), " %)")) %>% 
  arrange(compound_score) %>% 
  select(-compound_score, -ASV) %>% 
  knitr::kable(col.names = c('Genus', 
                             'Prev',
                             'Rank',
                             'Prev',
                             'Rank',
                             'ASV label'),
               align = 'lccccr') %>% 
  kableExtra::add_header_above(top_headers_names) %>% # Top header
  kableExtra::kable_styling(bootstrap_options = c('striped', 'hover'), 
                            full_width = FALSE) %>% 
  save_kable("2023/out/ASV_Genus_prevalence.html",
             self_contained = TRUE)

# Taxonomy table
joined_table %>% 
  select(ASV_label, ASV) %>% 
  left_join(melted %>% 
              select(OTU, Genus, Family, Order, Class, Phylum) %>% 
              unique,
            join_by(ASV == OTU)) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling('striped') %>% 
  save_kable("2023/out/ASV_labels.html",self_contained = TRUE)


