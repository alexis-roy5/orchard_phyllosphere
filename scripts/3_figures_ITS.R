## ITS figures ##

library(dplyr)
library(tidyverse)
library(phyloseq)
library(readxl)
library(rstatix)
library(ggpubr)
library(kableExtra)
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/myFunctions.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')


ps_ITS <- readRDS("2023/out/ps_ITS.rds")
ps.rarefied.ITS = rarefy_even_depth(ps_ITS, rngseed=1, sample.size=0.9*min(sample_sums(ps_ITS)), replace=F)

df.rarefied.ITS <- psmelt(ps.rarefied.ITS)

################
# ASV.SEQUENCE #
################

# Prep data
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

ps.stats <- ps.stats.prep(ps_ITS, "ITS") %>% 
  mutate(across(where(is.numeric), ~ format(round(., 0),big.mark=',')))

ps.stats.k <- kable(ps.stats, "html", align = "c") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(
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


html_file <- file.path("./2023/out/table_asv_sequence.html")
save_kable(ps.stats.k, file = html_file)

# Use webshot to convert the HTML file to PDF
webshot(html_file, 
        file.path("./2023/out/table_asv_sequence.html"),
        cliprect = "viewport")

###############
## BAR CHART ##
###############

library(RColorBrewer)
library(ggh4x) # for nested facets 

which_taxrank <- 'Order'

# Agglomerate ps object and melt into a long data frame
  # it cannot be in the function because we need
melted.sample <- ps.rarefied.ITS %>% 
  tax_glom(taxrank = which_taxrank) %>% 
  psmelt() %>% # "melts" the phyloseq object in a single table
  filter(Abundance != 0) %>% 
  group_by(cultivar, time, orchard, type, practice, !!sym(which_taxrank)) %>% # summing the abundance of replicates for weighted mean. 
  summarise(Abundance = sum(Abundance),
            .groups = 'drop') %>% 
  mutate(Sample = paste0(time,orchard,cultivar)) %>% # grouping per individuals of same time, cultivar and orchard.
  ungroup %>% 
  group_by(Sample) %>% # compute relative abundance by sample:
  mutate(relAb = Abundance/sum(Abundance))

# Find top taxa and create "Others" category
nTaxa <- 10
(top_taxa.sample <- topTaxa(melted.sample, which_taxrank, nTaxa)) # custom function

# Define nTaxa+2 levels (10 top taxa + Others and Unclassified categories)
# This will allow us to reorder factor levels
top_taxa_lvls.sample <- top_taxa.sample %>% 
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  arrange(relAb) %>% pull(aggTaxo) %>% 
  as.character %>% # Others first:
  setdiff(., c('Others', 'Unclassified')) %>% 
  c('Others', 'Unclassified', .)

# add a variable with only the top taxa identifiers
plot.df.sample <- melted.sample %>% # use a left-join :
  left_join(top_taxa.sample %>% select(-relAb), by = which_taxrank) %>%
  mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls.sample),
         time = factor(time, levels = c('May', 'July', 'August'))) # reorder factor levels

# Plot !
sample.chart <- plot.df.sample %>% 
  ggplot(aes(x = Sample, y = relAb, fill = aggTaxo)) +
  geom_col() +
  theme_light() +
  facet_nested(cols=vars(orchard,time), scales = 'free', space = 'free') +
  scale_fill_brewer(palette = 'Set3') +
  labs(fill = which_taxrank) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())

############################
# All samples across sites #
############################

site.chart <- plot.df.sample %>% 
  ggplot(aes(x = Sample, y = relAb, fill = aggTaxo)) +
  geom_col() +
  theme_light() +
  facet_nested(cols=vars(practice, orchard), scales = 'free', space = 'free') +
  scale_fill_brewer(palette = 'Set3') +
  labs(fill = which_taxrank) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank()); site.chart

########################### 
# All samples across time #
###########################

time.chart <- plot.df.sample %>% 
  ggplot(aes(x = Sample, y = relAb, fill = aggTaxo)) +
  geom_col() +
  theme_light() +
  facet_nested(cols=vars(time), scales = 'free', space = 'free') +
  scale_fill_brewer(palette = 'Set3') +
  labs(fill = which_taxrank) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank()) ; time.chart

###########################
# Conventional vs organic #
###########################

practice.chart <- plot.df.sample %>% 
  ggplot(aes(x = Sample, y = relAb, fill = aggTaxo)) +
  geom_col() +
  theme_light() +
  facet_nested(cols=vars(practice), scales = 'free', space = 'free') +
  scale_fill_brewer(palette = 'Set3') +
  labs(fill = which_taxrank) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())


################
## DIVERSITY ###
################



# PREPARING DATA
alpha.diversity.metrics.df.ITS <- df.rarefied.ITS %>%
  filter(Abundance>0) %>% # filtrer les lignes (ASV) de 0
  mutate(Sample = paste0(time,orchard,cultivar)) %>% # prendre en compte les réplicats
  group_by(Sample) %>% # grouper par sample
  mutate(Richness = n(), # unique ASV par sample
         RelativeAbundance = Abundance / sum(Abundance), 
         ShannonComponent = ifelse(RelativeAbundance > 0, 
                                   -RelativeAbundance * log(RelativeAbundance), # NOTE: fonctionnement du ifelse; ifelse(test, yes, no)
                                   0), # calculer pour chaque sample le (abondance relative * log(abondance relative)) terme
         SimpsonIndex = (sum(RelativeAbundance^2)), 
  ) %>%
  summarise( # keeping all my variables important
    code = dplyr::first(code),
    time = dplyr::first(time),
    orchard = dplyr::first(orchard),
    practice = dplyr::first(practice),
    type = dplyr::first(type),
    cultivar = dplyr::first(cultivar),
    Richness = dplyr::first(Richness), # Preserver mes variable avec dplyr::first(column_name)
    ShannonIndex = sum(ShannonComponent),  
    expShannonIndex = exp(sum(ShannonIndex)), # formule pour expShannonIndex
    SimpsonIndex = dplyr::first(SimpsonIndex),
    InverseSimpsonIndex = 1 / sum(RelativeAbundance^2), # formule pour InverseSimpsonIndex
  ) %>% 
  ungroup() 

# mettre en ordre inverse alphabétique la variable "time"
alpha.diversity.metrics.df.ITS <- alpha.diversity.metrics.df.ITS %>%
  mutate(time = factor(time, levels = c('May', 'July', 'August')))


# STATS #
t_test <- t_test(alpha.diversity.metrics.df.ITS, ShannonIndex ~ time)

# Vérification des conditions du test de t
# normalité des résidus; homoscédasticité des variances; indépendance des observation

# normalité
model.ITS <- lm(ShannonIndex ~ time, data = alpha.diversity.metrics.df.ITS) # ok?
residuals <- residuals(model.ITS)
hist(residuals) # visualisation
shapiro.test(residuals(model.ITS)) # shapiro.wilk


# variances
leveneTest(ShannonIndex ~ time, data = alpha.diversity.metrics.df.ITS)
# if p-value < 0.05 -> variances are not equal

# indépendance
plot(model.ITS$fitted.values, residuals)
abline(h = 0, col = "red")


my_comparisons <- list( c("May", "July"), c("July", "August"),c("May", "August") ) # for stats in my plots


# figures inspired from the manuscript 2023 16S of Sophie Boutin

plot.A <- ggplot(alpha.diversity.metrics.df.ITS, aes(x = time, y = ShannonIndex)) +
  geom_boxplot(aes(fill = practice), alpha = 0.50, outlier.shape = NA, size = 0.5) +  
  labs(
    x = NULL,
    y = "Shannon Index",
    subtitle = "Data from 2023",
    fill = "Practice"
  ) +
  facet_wrap(~code, ncol = 6) +
  #stat_compare_means(
    #comparisons = my_comparisons, 
    #label = "p.format", 
    #method = 'wilcox.test', 
    #p.adjust.method = "bonferroni") # Only show significant p-values
  #) + 
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  scale_fill_brewer(palette = "Set1"); plot.A

plot.B <- ggplot(alpha.diversity.metrics.df.ITS, aes(x = time, y = ShannonIndex)) +
  geom_boxplot(aes(fill = practice), alpha = 0.50, outlier.shape = NA, size = 0.5) +  
  labs(
    x = NULL,
    y = "Shannon Index",
    subtitle = "Data from 2023",
    fill = "Practice"
  ) +
  facet_wrap(~practice, ncol = 6) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", 
  method = 'wilcox.test', p.adjust.method = "bonferroni") + 
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()  
  ) +
  scale_fill_brewer(palette = "Set1")

plot.C <- ggplot(alpha.diversity.metrics.df.ITS, aes(x = time, y = ShannonIndex)) +
  geom_boxplot(aes(fill = practice), alpha = 0.50, outlier.shape = NA, size = 0.5) +  
  labs(
    x = NULL,
    y = "Shannon Index",
    subtitle = "Data from 2023",
    fill = "Practice"
  ) +
  # stat_compare_means(comparisons = my_comparisons, label = "p.format", 
  # method = 'wilcox.test', p.adjust.method = "bonferroni") + 
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Supprime le quadrillage principal
    panel.grid.minor = element_blank()  # Supprime le quadrillage mineur
  ) +
  scale_fill_brewer(palette = "Set1")

