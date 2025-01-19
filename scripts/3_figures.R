#######################
###### PHYLOSEQ #######
#######################


library(tidyverse)
library(phyloseq)
library(permute)
library(lattice)
library(vegan)
library(cowplot)

setwd("/Users/alexisroy/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere")

ps_16S

# raréfaction 
ps.rarefied = rarefy_even_depth(ps_16S, rngseed=1, sample.size=0.9*min(sample_sums(ps_16S)), replace=F)
# 0,9: pour réchantilloner le petit sampling. 

# faire fondre l'objet phyloseq
phylo.df <- psmelt(ps_16S) # pour comparaison
phylo.rarefied.df <- psmelt(ps.rarefied)

dim(phylo.df) # voir différences entre raréfier et non raréfier
dim(phylo.rarefied.df)

# estimer la richesse de tous les échantillons - nb ASVs
asv.nombre <- phylo.rarefied.df$OTU %>% 
  unique() %>% 
  length()
  # on aurait pu le faire en regardant seulement l'objet phyloseq

# estimer la richesse par sample
richesse.sample.df <- phylo.rarefied.df %>%
  filter(Abundance > 0) %>%             # Include only OTUs with abundance > 0
  group_by(Sample, time) %>%            # Group by `Sample` and `time`
  summarise(Richness = n_distinct(OTU), .groups = "drop") %>%  # Count unique OTUs
  left_join(phylo.rarefied.df %>% select(Sample, time, OTU, Abundance), 
            by = c("Sample", "time"))
  

# ABONDANCE totale
sample_sums <- phylo.rarefied.df %>%
  group_by(Sample) %>%              # Group by the Sample column
  summarize(TotalAbundance = sum(Abundance)) # Sum the Abundance column for each group

# ABONDANCE relative
phylo_with_rel_abundance <- richesse.sample.df %>%
  group_by(Sample) %>%                                    # Group by Sample
  mutate(TotalAbundance = sum(Abundance),                # Calculate total abundance for each sample
         RelativeAbundance = Abundance / TotalAbundance) %>% # Compute relative abundance
  ungroup() %>%  # Ensure ungrouping to avoid unintended grouping effects
  select(Sample, time, OTU, TotalAbundance, RelativeAbundance, Richness)

# Shannon index + exp(Shannon Index)
richness.shannon.df <- phylo_with_rel_abundance %>%
  group_by(Sample, time) %>%                           # Include 'time' in the grouping
  mutate(ShannonComponent = ifelse(RelativeAbundance > 0, 
                                   -RelativeAbundance * log(RelativeAbundance), 
                                   0)) %>%  # Calculate each p_i * log(p_i) term
  summarise(
    ShannonIndex = sum(ShannonComponent),     # Summarize Shannon components
    expShannonIndex = exp(sum(ShannonIndex)), 
    Richness = dplyr::first(Richness),    # Preserve pre-calculated Richness
    time = dplyr::first(time),            # Preserve 'time'
    RelativeAbundance = list(RelativeAbundance)  # Store Relative Abundance as a list
  ) %>%
  tidyr::unnest(cols = c(RelativeAbundance)) %>%  # Expand Relative Abundance
  ungroup()
# NOTE: fonctionnement du ifelse; ifelse(test, yes, no)

# inverse Gini-Simpson index
richness.shannon.simpson.df <- richness.shannon.df %>%
  group_by(Sample, time) %>%  # Group by Sample and time
  mutate(SimpsonIndex = sum(RelativeAbundance^2)) %>%  # Calculate SimpsonIndex
  summarise(
    Richness = max(Richness),  # selon moi, max fonctionne
    ShannonIndex = max(ShannonIndex),
    expShannonIndex =  dplyr::first(expShannonIndex),
    SimpsonIndex = max(SimpsonIndex),  # or use other appropriate aggregation
    InverseSimpsonIndex = 1 / sum(RelativeAbundance^2),  # Calculate the Inverse Simpson Index
    .groups = 'drop'  # Ungroup after summarizing
  ) %>%
  ungroup()

# mettre en ordre inverse alphabétique la variable "time"
richness.shannon.simpson.df <- richness.shannon.simpson.df %>%
  mutate(time = factor(time, levels = rev(sort(unique(time)))))

# violin plots

# modifer le plot - BRILLANT
R.S.S.index.long <- richness.shannon.simpson.df %>% 
  pivot_longer(cols = c(Richness, ShannonIndex, expShannonIndex, SimpsonIndex, InverseSimpsonIndex),
               names_to = "DiversityMeasure",
               values_to = "DiversityIndex") 
  

# filter for the normal diversity values
R.S.S.index.normal <- R.S.S.index.long %>%
  filter(DiversityMeasure %in% c("Richness", "ShannonIndex", "SimpsonIndex")) 

# filter for the Hill diversity values
R.S.S.index.hill <- R.S.S.index.long %>%
  filter(DiversityMeasure %in% c("Richness", "expShannonIndex", "InverseSimpsonIndex")) 


# mettre en ordre les mesures de diversité
R.S.S.index.normal$DiversityMeasure <- factor(R.S.S.index.normal$DiversityMeasure,
                                              levels = c("Richness", "ShannonIndex", "SimpsonIndex"))

R.S.S.index.hill$DiversityMeasure <- factor(R.S.S.index.hill$DiversityMeasure,
                                              levels = c("Richness", "expShannonIndex", "InverseSimpsonIndex"))




# time to plot 
  # normal values
plot.1 <- ggplot(R.S.S.index.normal, aes(x = time, y = DiversityIndex, color = time, fill = time)) +
  geom_boxplot(alpha = 0.25) +
  geom_jitter(width = 0.15, height = 0.05 , alpha = 0.5, size = 0.75) +
  labs(x = "Time", y = "Diversity Index",
       title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
       caption = "Data from 2023", color = "Time", fill = "Time") +
  facet_wrap(~ DiversityMeasure, scales = "free_y",
             labeller = as_labeller(c("Richness" = "Richness", 
                                      "ShannonIndex" = "Shannon Index", 
                                      "SimpsonIndex" = "Simpson Index"))) +
  theme_light() +
  theme(legend.position = "bottom")

# time to plot 
  # hill values
 plot.2 <- ggplot(R.S.S.index.hill, aes(x = time, y = DiversityIndex, color = time, fill = time)) +
  geom_boxplot(alpha = 0.25) +
  geom_jitter(width = 0.15, height = 0.05 , alpha = 0.5, size = 0.75) +
  labs(x = "Time", y = "Diversity Index",
       title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
       caption = "Data from 2023", color = "Time", fill = "Time") +
  facet_wrap(~ DiversityMeasure,
             labeller = as_labeller(c("Richness" = "Richness", 
                                      "expShannonIndex" = "Exponential Shannon Index", 
                                      "InverseSimpsonIndex" = "Inverse Simpson Index"))) +
  theme_light() +
  theme(legend.position = "bottom")

plot.1
plot.2
